//
// Created by luan on 10/09/17.
//

#ifndef STATIC_HEA_MINMIN_H
#define STATIC_HEA_MINMIN_H


#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>
#include "Data.h"
#include "Chromosome.h"
#include <boost/bimap.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/join.hpp>


using namespace TCLAP;
using namespace std;



// === MinMin Task scheduler === //

inline double ST(Data* data, int task, int vm, vector<double> ft_vector, vector<double> queue) {
    double max_pred_time = 0;
    for (auto tk : data->prec.find(task)->second)
        max_pred_time = std::max(max_pred_time, ft_vector[tk]);
    return std::max(max_pred_time, queue[vm]);
}

inline double transferTime(File file, VMachine vm1, VMachine vm2) {
    if (vm1.id != vm2.id) {
        auto link = std::min(vm1.bandwidth, vm2.bandwidth);
        return ceil(file.size / link);
    } else return 0;
}

inline double
FT(Data* data, Task task, VMachine vm, vector<double> &ft_vector, vector<double> &queue, vector<int> file_place, double lambda) {
    double start_time = 0;
    double read_time = 0;
    double write_time = 0;

    if (task.id != data->id_root && task.id != data->id_sink) {
        // Compute Start Time
        start_time = ST(data, task.id, vm.id, ft_vector, queue);
        // Read time;
        for (auto in : task.input) {

            auto file = data->file_map.find(in)->second;
            int vm_id = file.is_static ? file.static_vm : file_place[file.id];

            auto vmj = data->vm_map.find(file_place[vm_id])->second;
            read_time += ceil(transferTime(file, vm, vmj) + (file.size * lambda));
        }
        //write time
        for (auto out : task.output) {
            auto file = data->file_map.find(out)->second;
            write_time += ceil(file.size * (2 * lambda));
        }

    } else if (task.id == data->id_sink) {
        for (auto tk : data->prec.find(task.id)->second)
            start_time = std::max(start_time, ft_vector[tk]);
    }

    auto run_time = ceil(task.base_time * vm.slowdown);

    // cout << "task: " << data.idToString(task) << " vm: " << vm << endl;
    //cout << "start time: " << start_time << " read time: " << read_time << " Run time: " << run_time << " write_time: " << write_time << endl;

    return start_time + read_time + run_time + write_time;
}

void
schedule(Data* data, list<int> avail_tasks, vector<double> &ft_vector, vector<double> &queue, vector<int> &file_place,
         list<int> &task_ordering, double lambda) {

    double min_time;

    int min_vm = 0;
    vector<double> task_min_time(data->size, 0);
    vector<int> vm_min_time(data->size, 0);

    // while all task wasn't scheduled, do:
    while (!avail_tasks.empty()) {
        auto global_min_time = numeric_limits<double>::max();
        auto global_min_task = 0;
        // 1. Compute time phase
        for (auto task_id : avail_tasks) {
            // Compute the finish time off all tasks in each Vm
            min_time = numeric_limits<double>::max();

            auto task = data->task_map.find(task_id)->second;
            for (int j = 0; j < data->vm_size; j++) {
                auto vm = data->vm_map.find(j)->second;
                double time = FT(data, task, vm, ft_vector, queue, file_place, lambda);
                if (time < min_time) { // Get minimum time and minimum vm
                    min_time = time;
                    min_vm = j;
                }
            }
            if (global_min_time > min_time) {
                global_min_time = min_time;
                global_min_task = task.id;
            }
            task_min_time[task.id] = min_time;//Save the min_time of task
            vm_min_time[task.id] = min_vm;// and save the Vm with the min_time
        }


        auto r = vm_min_time[global_min_task];//r resource with min time in relation of min_task;
        //Update auxiliary structures (queue and ft_vector)
        ft_vector[global_min_task] = task_min_time[global_min_task];
        queue[r] = task_min_time[global_min_task];

        task_ordering.push_back(global_min_task);
        file_place[global_min_task] = r;

        //update file_place
        auto task = data->task_map.find(global_min_task)->second;
        for (auto file : task.output) {
            file_place[file] = r;
        }

        avail_tasks.remove(global_min_task);//remove task scheduled
    }

}


#endif //STATIC_HEA_MINMIN_H
