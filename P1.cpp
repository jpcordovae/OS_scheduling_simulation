#include <stdlib.h>
#include <string.h>
#include <chrono>
#include <thread>
#include <list>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <assert.h>
#include <iomanip>

#include <boost/tuple/tuple.hpp>
#include <boost/array.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>
#include <boost/bind.hpp>

#include "gnuplot-iostream.h"

#define ROUND_ROBIN_BITS 1
//#define P1_DEBUG

/// Global variables and structures
long time_bit = 1; // 1ms time slice
unsigned long long simulated_time = 0;

typedef struct st_process{
  pid_t pid;
  st_process *parent_process;
  st_process *next_process;
  st_process *back_process;
  st_process *child_process; // double linked list
  st_process *last_child_process; // last process of the child process
  unsigned int priority; // 2 to 16
  unsigned long long serviced_time; // time which has been run by the CPU
  unsigned long long service_time; // kernel process has 0 (zero) service time
  unsigned long long arrival_time; // time at it show up in the kernel process list
  unsigned long long remaining_time; // time left until finish work
  unsigned long long end_time; // timestamp when was taked out from running process
  unsigned long long wait_time; // guess what...
  unsigned long long turnarround_time; // how much time took the kernel to complete the process
  unsigned long long normalized_turnarround_time; // fos statistics
  int rr; // roud robin count

  st_process():parent_process(NULL),
               next_process(NULL), back_process(NULL),
               child_process(NULL),last_child_process(NULL) {}

  bool operator()(st_process *a, st_process *b) {
    return a->arrival_time < b->arrival_time;
  }
} st_process;

std::list<st_process*> process_list; // temporal conatiner with processes to add to kernel list
std::list<st_process*> finished_process_list;
st_process *priority_process[16];
/// functions
st_process *allocate_process(st_process *parent_process,long arrival_time, long max_service_time, pid_t &pid)
{
  st_process *tmp;
  tmp = (st_process*) malloc(sizeof(st_process));
  memset(tmp,0x00,sizeof(st_process));
  tmp->pid = pid++;
  tmp->priority = 1;
  tmp->arrival_time = (rand()%arrival_time)+1;
  tmp->service_time = (rand()%max_service_time)+1;
  tmp->serviced_time = 0L;
  tmp->wait_time = 0L;
  tmp->rr = 0;
  tmp->child_process = NULL;
  tmp->last_child_process = NULL;
  tmp->next_process = NULL;
  tmp->back_process = NULL;
  tmp->parent_process = NULL;
  return tmp;
}

/*st_process allocate_process(st_process *parent_process, unsigned long long arrival_time, unsigned long long max_service_time, pid_t &pid)
{
  st_process tmp;
  memset((void*)&tmp,0x00,sizeof(st_process));
  tmp.pid = pid++;
  tmp.priority = 1;
  tmp.rr = 0;
  tmp.arrival_time = (rand()%arrival_time)+1;
  tmp.service_time = (rand()%max_service_time)+1;
  tmp.parent_process = parent_process;
  return tmp;
  }*/

void delete_process_list(std::list<st_process*>&list_of_process)
{
  std::for_each(list_of_process.begin(),
                list_of_process.end(),
                [&](st_process *proc) {
                  st_process *tmp = proc; // secure delete
                  free(tmp);
                  proc = NULL;
                });
  list_of_process.clear();
}

st_process *add_process(st_process *parent_process, st_process *process)
{
  assert(process!=NULL);

  if(parent_process->child_process==NULL) {
    parent_process->child_process = process;
    parent_process->child_process->next_process = NULL;
    parent_process->child_process->back_process = NULL;
  } else {
    parent_process->last_child_process->next_process = process;
    process->back_process = parent_process->last_child_process;
  }
  parent_process->last_child_process = process;
  process->parent_process = parent_process;

  return process;
}

int count_child_process(st_process *p)
{
  int i = 0;
  if(!p) return 0;
  st_process *stP = p->child_process;
  if(!stP) return 0;
  i++;
  while(stP->next_process) {
    i++;
    stP = stP->next_process;
  }
  return i;
}

void swap_process(st_process *a, st_process *b)
{
  st_process *a_back = a->back_process;
  st_process *a_next = a->next_process;
  st_process *b_back = b->back_process;
  st_process *b_next = b->next_process;

  b->next_process = a_next;
  b->back_process = a_back;

  a->next_process = b_next;
  a->back_process = a_next;

  b_next->back_process = a;
  b_back->next_process = b;
}

// remove process from linked list and return, also fix the child_process pointer
st_process * remove_process(st_process *parent_process, st_process *process_to_delete)
{
  assert(process_to_delete);

  if(process_to_delete == parent_process->child_process) {
    parent_process->child_process = process_to_delete->next_process;
    process_to_delete->back_process = NULL;
    process_to_delete->next_process = NULL;
    return process_to_delete;
  }

  if(process_to_delete->next_process){
    process_to_delete->next_process->back_process = process_to_delete->back_process;
  }

  if(process_to_delete->back_process){
    process_to_delete->back_process->next_process = process_to_delete->next_process;
  }

  // fix last_child_process
  //if(process_to_delete==parent_process->last_child_process)
  //  parent_process->last_child_process = process_to_delete->back_process;

  //parent_process->last_child_process = parent_process->child_process;
  //while(parent_process->last_child_process != NULL) {
  //  parent_process->last_child_process = parent_process->last_child_process->next_process;
  //}

  process_to_delete->back_process = NULL;
  process_to_delete->next_process = NULL;

  return process_to_delete;
}

size_t count_process_with_priority(st_process* local_process_list,
                                   int priority)
{
  size_t counter = 0;
  st_process *tmp = local_process_list;
  while(tmp){
    if(tmp->priority==priority){
      counter++;
    }
    tmp = tmp->next_process;
  }
  return counter;
}

st_process **get_list_of_process_by_priority(st_process *local_process_list,
                                             int priority, int &counter)
{
  counter = count_process_with_priority(local_process_list,priority);
  if(!counter) return NULL;
  st_process **new_list = (st_process **) malloc(sizeof(st_process*)*counter);
  memset(new_list,0x00,sizeof(st_process*)*counter);
  st_process *tmp = local_process_list;
  int itmp=0;
  while(tmp!=NULL){
    if(tmp->priority==priority)
      new_list[itmp++] = tmp;
    tmp = tmp->next_process;
  }
  assert(counter==itmp);
  return new_list;
}

st_process *first_serve(st_process *parent_process, st_process *running_process)
{
  if(running_process!=NULL) {
    if(running_process->serviced_time>=running_process->service_time) {
      running_process->end_time = simulated_time;// save this information for analysis
      running_process->turnarround_time = running_process->service_time + running_process->wait_time;
      running_process->normalized_turnarround_time = running_process->turnarround_time / running_process->service_time;
      st_process *tmp = remove_process(parent_process,running_process);//remove the process froms the list and fix the child process of the parent if is the first in the list
#ifdef P1_DEBUG
      printf("proces %d with serviced time of %lld (%lld) removed at %lld (PC: %d).\n",
             tmp->pid, tmp->serviced_time, tmp->service_time, simulated_time,count_child_process(parent_process));
#endif
      finished_process_list.push_back(tmp); // save the process for future analysis
      running_process = NULL;
    }/*else {
      return running_process;
      }*/
  }

  for(int priority=1;priority<17;priority++){
    int number_of_process;
    st_process **p_list = get_list_of_process_by_priority(parent_process->child_process,priority,number_of_process);
    if(p_list) {
      running_process = p_list[0];
      free(p_list);
      return running_process;
    }
  }

  return running_process;
}

st_process *round_robin(st_process *parent_process, st_process *running_process)
{
  // check finish time
  if(running_process) {
    running_process->rr++;
    if(running_process->serviced_time == running_process->service_time){
      running_process->end_time = simulated_time;
      running_process->turnarround_time = running_process->service_time + running_process->wait_time;
      assert(running_process->serviced_time>0);
      running_process->normalized_turnarround_time = running_process->turnarround_time / running_process->service_time;
      st_process *tmp = remove_process(parent_process,running_process);
#ifdef P1_DEBUG
      printf("proces %d with serviced time of %lld (ST: %lld) removed at %lld (PC: %d).\n",
             tmp->pid, tmp->serviced_time, tmp->service_time,
             simulated_time,count_child_process(parent_process));
#endif
      finished_process_list.push_back(tmp);
    }
    // if running_process is null need to check next process by round robin index
    if(running_process->rr == ROUND_ROBIN_BITS) {
      running_process->rr = 0;
      // find next process
      // if running process has not finished, ned to check the round robin index
      for(int priority=1;priority<17;priority++){
        int number_of_process;
        st_process **p_list = get_list_of_process_by_priority(parent_process->child_process,
                                                              priority,
                                                              number_of_process);
        if(p_list) {
          for(int i=0;i<number_of_process;i++){
            if(p_list[i]==running_process){
              running_process = p_list[(i+1)%number_of_process];
              break;
            } //else {
              running_process = p_list[0];
              //}
          }
          free(p_list);
          return running_process;
        }
      }
    }

  } else {
    for(int priority=1;priority<17;priority++){
      int number_of_process;
      st_process **p_list = get_list_of_process_by_priority(parent_process->child_process,
                                                            priority,
                                                            number_of_process);
      if(p_list) {
        running_process = p_list[0];
        free(p_list);
      }
    }
  }

  return running_process;
}

st_process *shortest_remaining(st_process *parent_process, st_process *running_process)
{
  // check finish time
  if(running_process) {
    if(running_process->serviced_time == running_process->service_time){
      running_process->end_time = simulated_time;
      running_process->turnarround_time = running_process->service_time + running_process->wait_time;
      running_process->normalized_turnarround_time = running_process->turnarround_time / running_process->service_time;
      st_process *tmp=remove_process(parent_process,running_process);
#ifdef P1_DEBUG
      printf("proces %d with serviced time of %lld (%lld) removed at %lld (PC: %d).\n",
             tmp->pid, tmp->serviced_time, tmp->service_time, simulated_time,count_child_process(parent_process));
#endif
      finished_process_list.push_back(tmp);
    }
    running_process = NULL;
  }

  for(int priority=1;priority<17;priority++){
    int number_of_process;
    st_process **p_list = get_list_of_process_by_priority(parent_process->child_process,
                                                          priority,
                                                          number_of_process);
    if(p_list) {
      st_process *sr_tmp = p_list[0];
      if( sr_tmp && (number_of_process>1) ){
        for(int i=1;i<number_of_process-1;i++){
          if(p_list[i]->remaining_time < p_list[i+1]->remaining_time) {
            sr_tmp = p_list[i];
          }
        }
      }
      running_process = sr_tmp;
      free(p_list);
    }
  }
  return running_process;
}
void kernel_thread(st_process *kernel_process, st_process *(*schedule)(st_process *,st_process *))
{
  bool done = false;
  st_process *running_process=NULL;

  while(!done) {
    // update simulation
    simulated_time++;

    // update serviced time (CPU_time) for running process
    if(running_process) {
      running_process->serviced_time++;
    }

    // calculate wait time
    st_process *tmp = kernel_process->child_process;
    while(tmp!=NULL) {
      if(tmp!=running_process) {
        tmp->wait_time++;
      }
      tmp = tmp->next_process;
    }

    // calculate remaining time
    tmp = kernel_process->child_process;
    while(tmp) {
      //assert(tmp->serviced_time<=tmp->service_time);
      if(tmp->serviced_time>=tmp->service_time)
        tmp->remaining_time = 0;
      else
        tmp->remaining_time = tmp->service_time-tmp->serviced_time;
      tmp = tmp->next_process;
    }

    // process arrival from process list (future process) to kernel process list (active process)
    if(!process_list.empty()) { // the list is in order, thats why this way to deal with it
      // in case will be more than one process with same arrival time
      while(process_list.front()->arrival_time == simulated_time) {
        add_process(kernel_process,process_list.front());
        process_list.pop_front();
#ifdef P1_DEBUG
        st_process *tmpP = kernel_process->last_child_process;
        printf("new process arrival with PID: %d, arrival_time: %lld (RT: %lld), service time: %lld (pc: %d)\n",
               tmpP->pid,tmpP->arrival_time,simulated_time,tmpP->service_time,count_child_process(kernel_process));
#endif
        if(process_list.empty()) {
          break;
        }
      }
    }
    // selection of the running process
    running_process = schedule(kernel_process,running_process);
    // check end of the program
    if(process_list.empty() && kernel_process->child_process==NULL) done = true;
  }
  //printf("Finished \n");
}

bool compare(st_process *a, st_process *b)
{
  return a->arrival_time < b->arrival_time;
}

typedef struct st_process_statistics{
  enum type { FIRST_SERV, ROUND_ROBIN, MIN_EXEC_LEFT };
  double arrival_time_mean;
  double arrival_time_stdev;
  double turnarround_time_mean;
  double turnarround_time_stdev;
  double normalized_turnarround_mean;
  double normalized_turnarround_stdev;
  double service_time_mean;
  double service_time_stdev;
  double waiting_time_mean;
  double waiting_time_stdev;
  unsigned long long maximum_waiting_time;
  long max_arrival_time; // set of the simulation
  long max_process;
  long max_service_time;
  int bit_times;
} st_process_statistics;

std::ostream &operator<<(std::ostream &os, st_process_statistics &a)
{
  os << "execution of " << a.max_process << " process with service time of " << a.max_service_time << std::endl;
  os << "arrival mean: " << a.arrival_time_mean << std::endl;
  //  os << "arrival stdev: " << a.arrival_time_stdev << std::endl;
  os << "turnarround mean: " << a.turnarround_time_mean << std::endl;
  os << "normalized turnarround mean:" << a.normalized_turnarround_mean << std::endl;
  //os << "turnarround stdev: " << a.turnarround_time_stdev << std::endl;
  os << "service time mean: " << a.service_time_mean << std::endl;
  //os << "service ttime stdev: " << a.service_time_stdev << std::endl;
  os << "waiting time mean: " << a.waiting_time_mean << std::endl;
  os << "maximum waiting time: " << a.maximum_waiting_time << std::endl;
  os << std::endl;
  return os;
}

void calculate_mean_dev(std::vector<unsigned long long> &v_data, double &mean, double &stdev)
{
  mean =  std::accumulate(v_data.begin(),v_data.end(),0.0)/double(v_data.size());
  std::for_each(v_data.begin(), v_data.end(), [&](long tmp){
      stdev += (double(tmp)-mean)*(double(tmp)-mean);
    });
  stdev /= v_data.size()-1;
}

unsigned long long get_maximum_waiting_time(std::list<st_process*> &list_of_process)
{
  unsigned long long tmp = 0;
  std::for_each(list_of_process.begin(),
                list_of_process.end(),
                [&](st_process *proc){
                  if(proc->wait_time > tmp){
                    tmp = proc->wait_time;
                  }
                });
  return tmp;
}

st_process_statistics calculate_statistics(std::list<st_process*> &list_of_process)
{
  st_process_statistics st_ps;

  std::vector<unsigned long long> arrival_time;
  std::vector<unsigned long long> turnarround_time;
  std::vector<unsigned long long> normalized_turnarround;
  std::vector<unsigned long long> service_time;
  std::vector<unsigned long long> waiting_time;

  std::for_each(list_of_process.begin(),list_of_process.end(),[&](const st_process *proc){
      arrival_time.push_back(proc->arrival_time);
      turnarround_time.push_back(proc->turnarround_time);
      normalized_turnarround.push_back(proc->normalized_turnarround_time);
      service_time.push_back(proc->service_time);
      waiting_time.push_back(proc->wait_time);
    });

  calculate_mean_dev(arrival_time, st_ps.arrival_time_mean, st_ps.arrival_time_stdev);
  calculate_mean_dev(turnarround_time, st_ps.turnarround_time_mean, st_ps.turnarround_time_stdev);
  calculate_mean_dev(service_time, st_ps.service_time_mean, st_ps.service_time_stdev);
  calculate_mean_dev(waiting_time, st_ps.waiting_time_mean, st_ps.waiting_time_stdev);
  calculate_mean_dev(normalized_turnarround, st_ps.normalized_turnarround_mean, st_ps.normalized_turnarround_stdev);

  st_ps.maximum_waiting_time = get_maximum_waiting_time(list_of_process);

  return st_ps;
}

// reset a linked list to avoid allocation
void clean_process_list(std::list<st_process*>&plist)
{
  std::for_each(plist.begin(),
                plist.end(),
                [&](st_process *proc){
                  memset(proc,0x00,sizeof(st_process));
                });
}

void copy_process_list(std::list<st_process*> &dest, std::list<st_process*> &origin)
{
  if(!dest.empty()) delete_process_list(dest);
  std::for_each(origin.begin(),
                origin.end(),
                [&](st_process *proc){
                  st_process *tmp = (st_process*)malloc(sizeof(st_process));
                  memcpy(tmp,proc,sizeof(st_process));
                  dest.push_back(tmp);
                });
}

int main (void)
{
  std::vector<st_process_statistics> statistics_first_serve;
  std::vector<st_process_statistics> statistics_round_robin;
  std::vector<st_process_statistics> statistics_shortest;

  std::list<st_process*> p_list_problem_1;
  std::list<st_process*> p_list_problem_2;
  std::list<st_process*> p_list_problem_3;

  st_process *kernel_process = (st_process*) malloc(sizeof(st_process));
  long max_process = 1000;
  pid_t pid_index = 0;
  long max_service_time = 100;
  long max_arrival_time = 100;

  // PROBLEM 1
  std::cout << std::endl << "PROBLEM 1" << std::endl;
  //  st_process *(*fp[3])(st_process *,st_process *, std::list<st_process*>&, std::list<st_process*>&) = { first_serve, round_robin, shortest_remaining };
  // kt[0] = first_serve;
  // kt[1] = round_robin;
  // kt[2] = shortest_remaining;
  std::cout << "alocating 1000 process (process list)" << std::endl;
  std::cout << "max arrival time: " << max_arrival_time << std::endl;
  std::cout << "max service time: " << max_service_time << std::endl;
  pid_index=0;
  for(int i=0;i<max_process; i++) {
    p_list_problem_1.push_back(allocate_process(kernel_process,max_arrival_time,max_service_time,pid_index));
  }

  copy_process_list(p_list_problem_2, p_list_problem_1);
  copy_process_list(p_list_problem_3, p_list_problem_1);

  for(size_t i=10;i<=100;i+=10) {
    //std::list<st_process*> process_list;
    //std::list<st_process*> finished_process_list;
    // first serve scheduling
    std::cout << "--------------------------------------------------" << std::endl;
    memset(kernel_process,0x00,sizeof(st_process));
    max_service_time = i;
    simulated_time = 0; // global variable
    delete_process_list(process_list);
    delete_process_list(finished_process_list);
    copy_process_list(process_list,p_list_problem_1);
    std::cout << "udpating service time between 1 and " << max_service_time << std::endl;
    for(st_process *proc : process_list){
      proc->service_time = rand()%(max_service_time-1)+1;
    }
    copy_process_list(p_list_problem_1, process_list); // backing up the process list for next shcedules types
    std::cout << "sorting proces list" << std::endl;
    process_list.sort(compare);
    std::cout << "simulation data:" <<
      std::endl << "   max_service_time: " << max_service_time <<
      std::endl << "   max_arrival_time: " << max_arrival_time <<
      std::endl << "   max process:      " << max_process      <<
      std::endl;
    std::cout << "runnig thread for first serve scheduling" << std::endl;
    std::thread kThread_hw(kernel_thread,kernel_process,first_serve);
    kThread_hw.join();
    std::cout << "calculating statistics" << std::endl;
    statistics_first_serve.push_back(calculate_statistics(finished_process_list));
    statistics_first_serve[statistics_first_serve.size()-1].max_process = max_process;
    statistics_first_serve[statistics_first_serve.size()-1].max_service_time = i;
    statistics_first_serve[statistics_first_serve.size()-1].max_arrival_time = max_arrival_time;

    // round robin scheduling
    memset(kernel_process,0x00,sizeof(st_process));
    simulated_time = 0; // global variable
    //pid_index = 0;
    delete_process_list(finished_process_list);
    delete_process_list(process_list);
    std::cout << "restoring process list" << std::endl;
    copy_process_list(process_list, p_list_problem_1);
    std::cout << "sorting process list" << std::endl;
    process_list.sort(compare);
    std::cout << "runnig thread for round robin scheduling" << std::endl;
    std::thread kThread_rr(kernel_thread,kernel_process,round_robin);
    kThread_rr.join();
    std::cout << "calculating statistics" << std::endl;
    statistics_round_robin.push_back(calculate_statistics(finished_process_list));
    statistics_round_robin[statistics_round_robin.size()-1].max_process = max_process;
    statistics_round_robin[statistics_round_robin.size()-1].max_service_time = i;
    statistics_round_robin[statistics_round_robin.size()-1].max_arrival_time = max_arrival_time;
    // shortest-remaining scheduling
    memset(kernel_process,0x00,sizeof(st_process));
    simulated_time = 0; // global variable
    //pid_index = 0;
    delete_process_list(finished_process_list);
    delete_process_list(process_list);
    std::cout << "restoring process list" << std::endl;
    copy_process_list(process_list, p_list_problem_1);
    std::cout << "sorting process list" << std::endl;
    process_list.sort(compare);
    std::cout << "runnig thread for shortest remaining scheduling" << std::endl;
    std::thread kThread_sh(kernel_thread,kernel_process,shortest_remaining);
    kThread_sh.join();
    std::cout << "calculating statistics" << std::endl;
    statistics_shortest.push_back(calculate_statistics(finished_process_list));
    statistics_shortest[statistics_shortest.size()-1].max_process = max_process;
    statistics_shortest[statistics_shortest.size()-1].max_service_time = i;
    statistics_shortest[statistics_shortest.size()-1].max_arrival_time = max_arrival_time;
  }

  std::cout << std::endl << "STATISTICS RESUME PROBLEM 1" << std::endl;
  std::cout << "FIRST SERVE" << std::endl;
    std::for_each(statistics_first_serve.begin(),
    statistics_first_serve.end(),
    [&](st_process_statistics st){
    std::cout << st << std::endl;
    });

    std::cout << "ROUND ROBIN" << std::endl;
    std::for_each(statistics_round_robin.begin(),
    statistics_round_robin.end(),
    [&](st_process_statistics st){
    std::cout << st << std::endl;
    });

    std::cout << "SHORTEST" << std::endl;
    std::for_each(statistics_shortest.begin(),
    statistics_shortest.end(),
    [&](st_process_statistics st){
    std::cout << st << std::endl;
    });

  // table output problem 1
  std::cout << std::setw(17) << " "
            << std::setw(30) << "arrival time"
            << std::endl;

  std::cout << std::setw(17) << "service time (ms)"
            << std::setw(10) << "FCFS" << " | "
            << std::setw(10) << "RR" << " | "
            << std::setw(10) << "SRT" << std::endl;

  for(size_t i=0;i<statistics_first_serve.size();i++) {
    std::cout << std::setw(17) << 10+i*10 <<
      std::setw(10) << statistics_first_serve[i].arrival_time_mean << " | " <<
      std::setw(10) << statistics_round_robin[i].arrival_time_mean << " | " <<
      std::setw(10) << statistics_shortest[i].arrival_time_mean << std::endl;
  }

  std::cout << std::endl;

  std::cout << std::setw(17) << " "
            << std::setw(30) << "normalized turnarround mean"
            << std::endl;

  std::cout << std::setw(17) << "service time (ms)"
            << std::setw(10) << "FCFS" << " | "
            << std::setw(10) << "RR" << " | "
            << std::setw(10) << "SRT"
            << std::endl;

  for(size_t i=0;i<statistics_first_serve.size();i++) {
    std::cout << std::setw(17) << 10+i*10 <<
      std::setw(10) << statistics_first_serve[i].normalized_turnarround_mean << " | " <<
      std::setw(10) << statistics_round_robin[i].normalized_turnarround_mean << " | " <<
      std::setw(10) << statistics_shortest[i].normalized_turnarround_mean <<
      std::endl;
  }

  std::cout << std::endl;

  std::cout << std::setw(17) << " "
            << std::setw(30) << "maximum wait time" << std::endl;

  std::cout << std::setw(17) << "service time (ms)"
            << std::setw(10) << "FCFS" << " | "
            << std::setw(10) << "RR" << " | "
            << std::setw(10) << "SRT" << std::endl;

  for(size_t i=0;i<statistics_first_serve.size();i++) {
    std::cout << std::setw(17) << 10+i*10 <<
      std::setw(10) << statistics_first_serve[i].maximum_waiting_time << " | " <<
      std::setw(10) << statistics_round_robin[i].maximum_waiting_time << " | " <<
      std::setw(10) << statistics_shortest[i].maximum_waiting_time <<
      std::endl;
  }

  //Gnuplot gp("gplot_script.gp");
  //gp << ""
  ////////////////////////////////////////////////////////////////////////////////////
  //
  // PROBLEM 2
  //
  ////////////////////////////////////////////////////////////////////////////////////
  std::cout << std::endl << "PROBLEM 2" << std::endl;
  std::cout << "---------------------------------------------------------" << std::endl;
  std::cout << "setting service_time: 80-100 (50%) and 1-10 (50%)" << std::endl;
  //clean_process_list(p_list_problem_2);
  size_t mid_list = std::floor(p_list_problem_2.size()/2);
  size_t for_index = 0;
  std::for_each(p_list_problem_2.begin(),
                p_list_problem_2.end(),
                [&](st_process *proc){
                  if(for_index<=mid_list)
                    proc->service_time = 80 + rand()%21;
                  else{
                    proc->service_time = 1 + rand()%10;
                  }
                  for_index++;
                });
  std::cout << "deleting process list" << std::endl;
  delete_process_list(process_list);
  std::cout << "deleting finished process list" << std::endl;
  delete_process_list(finished_process_list);
  std::cout << "restoring process list" << std::endl;
  copy_process_list(process_list,p_list_problem_2);
  std::cout << "cleaning statistics containers" << std::endl;
  statistics_first_serve.clear();
  statistics_round_robin.clear();
  statistics_shortest.clear();
  // first serve scheduling
  memset(kernel_process,0x00,sizeof(st_process));
  simulated_time = 0;
  //pid_index = 0;
  std::cout << "sorting proces list" << std::endl;
  process_list.sort(compare);
  std::cout << "simulation data:" <<
    std::endl << "service_time: 80-100 (50%) and 1-10 (50%)"  <<
    std::endl << "max_arrival_time: " << max_arrival_time <<
    std::endl;
  std::cout << "runnig thread for first serve scheduling" << std::endl;
  std::thread kThread_p2_fs(kernel_thread,kernel_process,first_serve);
  kThread_p2_fs.join();
  std::cout << "calculating statistics" << std::endl;
  statistics_first_serve.push_back(calculate_statistics(finished_process_list));
  statistics_first_serve[statistics_first_serve.size()-1].max_process = max_process;
  statistics_first_serve[statistics_first_serve.size()-1].max_service_time = 0;
  statistics_first_serve[statistics_first_serve.size()-1].max_arrival_time = max_arrival_time;

  std::list<st_process*> fs_low_priority_servicetime;
  std::list<st_process*> fs_high_priority_servicetime;
  for(st_process *proc: finished_process_list) {
    if(proc->service_time<=10) fs_high_priority_servicetime.push_back(proc);
    else fs_low_priority_servicetime.push_back(proc);
  }
  st_process_statistics fs_low_priority_statistics = calculate_statistics(fs_low_priority_servicetime);
  st_process_statistics fs_high_priority_statistics = calculate_statistics(fs_high_priority_servicetime);
  // round robin scheduling
  std::cout << "---------------------------------------------------------" << std::endl;
  memset(kernel_process,0x00,sizeof(st_process));
  simulated_time = 0; // global variable
  std::cout << "deleting process list" << std::endl;
  delete_process_list(process_list);
  std::cout << "deleting finished process list" << std::endl;
  delete_process_list(finished_process_list);
  std::cout << "restoring process list" << std::endl;
  copy_process_list(process_list,p_list_problem_2);
  std::cout << "sorting proces list" << std::endl;
  process_list.sort(compare);
  std::cout << "runnig thread for round robin scheduling" << std::endl;
  std::thread kThread_p2_rr(kernel_thread,kernel_process,round_robin);
  kThread_p2_rr.join();
  std::cout << "calculating statistics" << std::endl;
  statistics_round_robin.push_back(calculate_statistics(finished_process_list));
  statistics_round_robin[statistics_round_robin.size()-1].max_process = max_process;
  statistics_round_robin[statistics_round_robin.size()-1].max_service_time = 0;
  statistics_round_robin[statistics_round_robin.size()-1].max_arrival_time = max_arrival_time;

  std::list<st_process*> rr_low_priority_servicetime;
  std::list<st_process*> rr_high_priority_servicetime;
  for(st_process *proc: finished_process_list) {
    if(proc->service_time<=10) rr_high_priority_servicetime.push_back(proc);
    else rr_low_priority_servicetime.push_back(proc);
  }
  st_process_statistics rr_low_priority_statistics = calculate_statistics(rr_low_priority_servicetime);
  st_process_statistics rr_high_priority_statistics = calculate_statistics(rr_high_priority_servicetime);
  // shortest-remaining scheduling
  std::cout << "---------------------------------------------------------" << std::endl;
  memset(kernel_process,0x00,sizeof(st_process));
  simulated_time = 0; // global variable
  std::cout << "deleting process list" << std::endl;
  delete_process_list(process_list);
  std::cout << "deleting finished process list" << std::endl;
  delete_process_list(finished_process_list);
  std::cout << "restoring process list" << std::endl;
  copy_process_list(process_list,p_list_problem_2);
  std::cout << "sorting proces list" << std::endl;
  process_list.sort(compare);
  std::cout << "runnig thread for shortest remaining scheduling" << std::endl;
  std::thread kThread_p2_sh(kernel_thread,kernel_process,shortest_remaining);
  kThread_p2_sh.join();
  std::cout << "calculating statistics" << std::endl;
  statistics_shortest.push_back(calculate_statistics(finished_process_list));
  statistics_shortest[statistics_shortest.size()-1].max_process = max_process;
  statistics_shortest[statistics_shortest.size()-1].max_service_time = 0;
  statistics_shortest[statistics_shortest.size()-1].max_arrival_time = max_arrival_time;
  std::list<st_process*> sr_low_priority_servicetime;
  std::list<st_process*> sr_high_priority_servicetime;
  for(st_process *proc: finished_process_list) {
    if(proc->service_time<=10) sr_high_priority_servicetime.push_back(proc);
    else sr_low_priority_servicetime.push_back(proc);
  }
  st_process_statistics sr_low_priority_statistics = calculate_statistics(sr_low_priority_servicetime);
  st_process_statistics sr_high_priority_statistics = calculate_statistics(sr_high_priority_servicetime);

  // table output
  std::cout << "turnarround time mean overall " << std::endl;
  std::cout << std::setw(17) << "service time (ms)"
            << std::setw(10) << "FCFS" << " | "
            << std::setw(10) << "RR" <<   " | "
            << std::setw(10) << "SRT" << std::endl;
  std::cout << std::setw(17) << "1ms - 100ms"
            << std::setw(10) << statistics_first_serve[0].turnarround_time_mean << " | "
            << std::setw(10) << statistics_round_robin[0].turnarround_time_mean << " | "
            << std::setw(10) << statistics_shortest[0].turnarround_time_mean << std::endl;

  std::cout << "turnarround time mean low service time " << std::endl;
  std::cout << std::setw(17) << "service time (ms)"
            << std::setw(10) << "FCFS" << " | "
            << std::setw(10) << "RR" <<   " | "
            << std::setw(10) << "SRT" << std::endl;
  std::cout << std::setw(17) << "1ms - 10ms"
            << std::setw(10) << fs_low_priority_statistics.turnarround_time_mean << " | "
            << std::setw(10) << rr_low_priority_statistics.turnarround_time_mean << " | "
            << std::setw(10) << sr_low_priority_statistics.turnarround_time_mean << std::endl;

  std::cout << std::endl;
  std::cout << "turnarround time mean high service time " << std::endl;
  std::cout << std::setw(17) << "service time (ms)"
            << std::setw(10) << "FCFS" << " | "
            << std::setw(10) << "RR" << " | "
            << std::setw(10) << "SRT" << std::endl;
  std::cout << std::setw(17) << "80ms - 100ms"
            << std::setw(10) << fs_high_priority_statistics.turnarround_time_mean << " | "
            << std::setw(10) << rr_high_priority_statistics.turnarround_time_mean << " | "
            << std::setw(10) << sr_high_priority_statistics.turnarround_time_mean << std::endl;

  // split by service time
  /*std::list<st_process*>lst_serv_time_process; // low service time
  std::list<st_process*>hst_serv_time_process; // high service time

  st_process_statistics st_s_lst;
  st_process_statistics st_s_hst;

  for(st_process *proc: finished_process_list) {
    if(proc->service_time<= 10) {
      lst_serv_time_process.push_back(proc);
    } else {
      hst_serv_time_process.push_back(proc);
    }
  }

  st_s_lst = calculate_statistics(lst_serv_time_process);
  st_s_hst = calculate_statistics(hst_serv_time_process);

  std::cout << std::endl;
  std::cout << "mean turnaround for each group" << std::endl;
  std::cout << "------------------------------" << std::endl;
  std::cout << std::setw(17) << " 1ms - 10ms"
            << std::setw(10) << st_s_lst.normalized_turnarround_mean << " | "
            << std::endl;
  std::cout << std::setw(17) << "80ms - 100ms"
            << std::setw(10) << st_s_lst.normalized_turnarround_mean << " | "
            << std::endl;*/

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // PROBLEM 3
  //
  ////////////////////////////////////////////////////////////////////////////////////
  std::cout << std::endl << "PROBLEM 3" << std::endl;
  std::cout << "---------------------------------------------------------" << std::endl;
  std::cout << "cleaning kernel" << std::endl;
  memset(kernel_process,0x00,sizeof(st_process));
  std::cout << "setting service_time: 80-100 (50%) and 1-10 (50%)" << std::endl;
  //clean_process_list(p_list_problem_2);
  mid_list = max_process*0.9;
  for_index = 0;
  copy_process_list(p_list_problem_3,p_list_problem_2);
  std::list<st_process*> low_priority;
  std::list<st_process*> high_priority;
  std::cout << "cleaning statistics containers" << std::endl;
  statistics_first_serve.clear();
  statistics_round_robin.clear();
  statistics_shortest.clear();
  std::cout << "fixing the priority of the list" <<  std::endl;
  std::for_each(p_list_problem_3.begin(),
                p_list_problem_3.end(),
                [&](st_process *proc){
                  if(for_index < mid_list) {
                    proc->priority = 1;
                    high_priority.push_back(proc);
                  } else{
                    proc->priority = 2+rand()%15;
                    low_priority.push_back(proc);
                  }
                  for_index++;
                });
  // first serve scheduling
  memset(kernel_process,0x00,sizeof(st_process));
  simulated_time = 0;
  std::cout << "deleting process list" << std::endl;
  delete_process_list(process_list);
  std::cout << "deleting finished process list" << std::endl;
  delete_process_list(finished_process_list);
  std::cout << "restoring process list" << std::endl;
  copy_process_list(process_list,p_list_problem_3);
  std::cout << "sorting proces list" << std::endl;
  process_list.sort(compare);
  std::cout << "simulation data:" << std::endl
            << "  service_time: 80-100 (50%) and 1-10 (50%)"  << std::endl
            << "  priority 1 (90%) and 2-16 (10%)" << std ::endl
            << "  max_arrival_time: " << max_arrival_time << std::endl;
  std::cout << "runnig thread for first serve scheduling" << std::endl;
  std::thread kThread_p3_fs(kernel_thread,kernel_process,first_serve);
  kThread_p3_fs.join();
  std::cout << "calculating statistics" << std::endl;
  statistics_first_serve.push_back(calculate_statistics(finished_process_list));
  statistics_first_serve[statistics_first_serve.size()-1].max_process = max_process;
  statistics_first_serve[statistics_first_serve.size()-1].max_service_time = 0;
  statistics_first_serve[statistics_first_serve.size()-1].max_arrival_time = max_arrival_time;

  std::list<st_process*> hp_list_fs; // high process list
  std::list<st_process*> lp_list_fs; // low priority list
  for(st_process *proc: finished_process_list) {
    if(proc->priority==1) hp_list_fs.push_back(proc);
    else lp_list_fs.push_back(proc);
  }
  st_process_statistics hp_stats_fs = calculate_statistics(hp_list_fs);
  st_process_statistics lp_stats_fs = calculate_statistics(lp_list_fs);
  // round robin scheduling
  std::cout << "---------------------------------------------------------" << std::endl;
  memset(kernel_process,0x00,sizeof(st_process));
  simulated_time = 0; // global variable
  std::cout << "deleting process list" << std::endl;
  delete_process_list(process_list);
  std::cout << "deleting finished process list" << std::endl;
  delete_process_list(finished_process_list);
  std::cout << "restoring process list" << std::endl;
  copy_process_list(process_list,p_list_problem_3);
  std::cout << "sorting proces list" << std::endl;
  process_list.sort(compare);
  std::cout << "simulation data:" << std::endl
            << "  service_time: 80-100 (50%) and 1-10 (50%)"  << std::endl
            << "  priority 1 (90%) and 2-16 (10%)" << std ::endl
            << "  max_arrival_time: " << max_arrival_time << std::endl;
  std::cout << "runnig thread for round robin scheduling" << std::endl;
  std::thread kThread_p3_rr(kernel_thread,kernel_process,round_robin);
  kThread_p3_rr.join();
  std::cout << "calculating statistics" << std::endl;
  statistics_round_robin.push_back(calculate_statistics(finished_process_list));
  statistics_round_robin[statistics_round_robin.size()-1].max_process = max_process;
  statistics_round_robin[statistics_round_robin.size()-1].max_service_time = 0;
  statistics_round_robin[statistics_round_robin.size()-1].max_arrival_time = max_arrival_time;

  std::list<st_process*> hp_list_rr; // high process list
  std::list<st_process*> lp_list_rr; // low priority list
  for(st_process *proc: finished_process_list) {
    if(proc->priority<=10) hp_list_rr.push_back(proc);
    else lp_list_rr.push_back(proc);
  }
  st_process_statistics hp_stats_rr = calculate_statistics(hp_list_rr);
  st_process_statistics lp_stats_rr = calculate_statistics(lp_list_rr);
// shortest-remaining scheduling
  std::cout << "---------------------------------------------------------" << std::endl;
  memset(kernel_process,0x00,sizeof(st_process));
  simulated_time = 0; // global variable
  std::cout << "deleting process list" << std::endl;
  delete_process_list(process_list);
  std::cout << "deleting finished process list" << std::endl;
  delete_process_list(finished_process_list);
  std::cout << "restoring process list" << std::endl;
  copy_process_list(process_list,p_list_problem_3);
  std::cout << "sorting proces list" << std::endl;
  process_list.sort(compare);
  std::cout << "simulation data:" << std::endl
            << "  service_time: 80-100 (50%) and 1-10 (50%)"  << std::endl
            << "  priority 1 (90%) and 2-16 (10%)" << std ::endl
            << "  max_arrival_time: " << max_arrival_time << std::endl;
  std::cout << "runnig thread for shortest remaining scheduling" << std::endl;
  std::thread kThread_p3_sh(kernel_thread,kernel_process,shortest_remaining);
  kThread_p3_sh.join();
  std::cout << "calculating statistics" << std::endl;
  statistics_shortest.push_back(calculate_statistics(finished_process_list));
  statistics_shortest[statistics_shortest.size()-1].max_process = max_process;
  statistics_shortest[statistics_shortest.size()-1].max_service_time = 0;
  statistics_shortest[statistics_shortest.size()-1].max_arrival_time = max_arrival_time;

  std::list<st_process*> hp_list_sr; // high process list
  std::list<st_process*> lp_list_sr; // low priority list
  for(st_process *proc: finished_process_list) {
    if(proc->priority==1) hp_list_sr.push_back(proc);
    else lp_list_sr.push_back(proc);
  }
  st_process_statistics hp_stats_sr = calculate_statistics(hp_list_sr);
  st_process_statistics lp_stats_sr = calculate_statistics(lp_list_sr);
  // table output
  std::cout << std::setw(17) << " " << std::setw(30) << "mean turnarround time" << std::endl;
  std::cout << std::setw(17) << "service time (ms)" <<
    std::setw(10) << "FCFS" << " | " <<
    std::setw(10) << "RR" <<   " | " <<
    std::setw(10) << "SRT" << std::endl;

  for(size_t i=0;i<statistics_first_serve.size();i++) {
    std::cout << std::setw(17) << 10+i*10 <<
      std::setw(10) << statistics_first_serve[i].turnarround_time_mean << " | " <<
      std::setw(10) << statistics_round_robin[i].turnarround_time_mean << " | " <<
      std::setw(10) << statistics_shortest[i].turnarround_time_mean << std::endl;
  }

  std::cout << std::endl;

  std::cout << std::setw(17) << " " << std::setw(30) << "normalized turnarround mean" << std::endl;
  std::cout << std::setw(17) << "service time (ms)"
            << std::setw(10) << "FCFS" << " | "
            << std::setw(10) << "RR" << " | "
            << std::setw(10) << "SRT" << std::endl;

  for(size_t i=0;i<statistics_first_serve.size();i++) {
    std::cout << std::setw(17) << 10+i*10 <<
      std::setw(10) << statistics_first_serve[i].normalized_turnarround_mean << " | " <<
      std::setw(10) << statistics_round_robin[i].normalized_turnarround_mean << " | " <<
      std::setw(10) << statistics_shortest[i].normalized_turnarround_mean <<
      std::endl;
  }
  //---------------------------
  std::cout << std::endl;
  std::cout << "statistic resume high priority process (value of 1)" << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << std::setw(23) << "ITEM"
            << std::setw(10) << "FCFS" << " | "
            << std::setw(10) << "RR"   << " | "
            << std::setw(10) << "SRT" << std::endl;
  std::cout << std::setw(23) << "normalized turnarround"
            << std::setw(10) << hp_stats_fs.normalized_turnarround_mean << " | "
            << std::setw(10) << hp_stats_rr.normalized_turnarround_mean << " | "
            << std::setw(10) << hp_stats_sr.normalized_turnarround_mean << std::endl;
  std::cout << std::setw(23) << "waiting time mean"
            << std::setw(10) << hp_stats_fs.waiting_time_mean << " | "
            << std::setw(10) << hp_stats_rr.waiting_time_mean << " | "
            << std::setw(10) << hp_stats_sr.waiting_time_mean << std::endl;
  //---------------------------
  std::cout << std::endl;
  std::cout << "statistic resume low priority process (values 2-16)" << std::endl;
  std::cout << "--------------------------------------------------" << std::endl;
  std::cout << std::setw(23) << "ITEM"
            << std::setw(10) << "FCFS" << " | "
            << std::setw(10) << "RR"   << " | "
            << std::setw(10) << "SRT" << std::endl;
  std::cout << std::setw(23) << "normalized turnarround"
            << std::setw(10) << lp_stats_fs.normalized_turnarround_mean << " | "
            << std::setw(10) << lp_stats_rr.normalized_turnarround_mean << " | "
            << std::setw(10) << lp_stats_sr.normalized_turnarround_mean << std::endl;
  std::cout << std::setw(23) << "waiting time mean"
            << std::setw(10) << lp_stats_fs.waiting_time_mean << " | "
            << std::setw(10) << lp_stats_rr.waiting_time_mean << " | "
            << std::setw(10) << lp_stats_sr.waiting_time_mean << std::endl;

}
