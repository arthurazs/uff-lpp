#ifndef HEA_HEFT_H
#define HEA_HEFT_H

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


// === HEFT ===//

/*Event represents the start and end time of a task(id)*/
struct Event {
    int id;
    double start = 0;
    double end = 0;
};


typedef map<int, vector<Event>> event_map;

vector<int> instersection(vector<int> v1, vector<int> v2) {
    vector<int> v3;

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));

    return v3;
}

// get makespan
double makespan(event_map orders) {
    auto mksp = 0.0;
    for (auto it : orders) {
        if (it.second.size() > 0)
            mksp = std::max(it.second.back().end, mksp);
    }
    return mksp;
}

double commcost_static(int id_task, int id_vm, Data* data, double lambda) {
    double cost = 0.0;

    auto task = data->task_map.find(id_task)->second;
    auto vm = data->vm_map.find(id_vm)->second;

    for (auto id_file : task.input) {
        auto file = data->file_map.find(id_file)->second;
        if (file.is_static && file.static_vm != vm.id) {
            auto s_vm = data->vm_map.find(file.static_vm)->second;
            auto bandwidth = std::min(s_vm.bandwidth, vm.bandwidth);
            cost += ceil((file.size / bandwidth) + (file.size * lambda));
        }
    }

    for (auto id_file : task.output) {//for each file write by task, do
        auto file = data->file_map.find(id_file)->second;
        cost += ceil(file.size * (2 * lambda));
    }

    return cost;
}

/*Compute communication cost of dynamic files*/
double commcost_dynamic(int id_taski, int id_taskj, int id_vmi, int id_vmj, Data* data, double lambda) {
    double cost = 0.0;

    if (id_vmi == id_vmj)
        return cost;

    auto task_i = data->task_map.find(id_taski)->second;
    auto task_j = data->task_map.find(id_taskj)->second;

    auto vm_i = data->vm_map.find(id_vmi)->second;
    auto vm_j = data->vm_map.find(id_vmj)->second;

    //get the lowest bandwidth
    auto bandwidth = std::min(vm_i.bandwidth, vm_j.bandwidth);

    //get files write by task_i and read by task_j
    auto vet_files = instersection(task_i.output, task_j.input);

    for (auto id_file : vet_files) {//for each file reading by task_j, do
        auto file = data->file_map.find(id_file)->second;//get file
        //if file is write by task_i and read by task_j, do
        cost += ceil((file.size / bandwidth) + (file.size * lambda));
    }

    return cost;
}

double compcost(int id_task, int id_vm, Data* data) {
    auto task = data->task_map.find(id_task)->second;
    auto vm = data->vm_map.find(id_vm)->second;
    return ceil(task.base_time * vm.slowdown);

}

/*average computation cost*/
double wbar(int id_task, Data* data) {
    double wbar_cost = 0.0;
    for (auto it : data->vm_map)
        wbar_cost += compcost(id_task, it.first, data);
    return wbar_cost / double(data->vm_size);
}

/*average communication cost*/
double cbar(int id_taski, int id_taskj, Data* data, double lambda) {
    double cbar_cost = 0.0;
    if (data->vm_size == 1)
        return cbar_cost;
    //get number of pairs
    auto n_pairs = data->vm_size * (data->vm_size - 1);

    //for each vm1, compute average static file communication cost
    for (auto vm1 : data->vm_map)
        cbar_cost = commcost_static(id_taskj, vm1.first, data, lambda);

    //for each pair of vms compute the average communication between taski and taskj
    for (auto vm1 : data->vm_map) {
        for (auto vm2 : data->vm_map) {
            if (vm1.first != vm2.first)
                cbar_cost += commcost_dynamic(id_taski, id_taskj, vm1.first, vm2.first, data, lambda);
        }
    }
    return 1. * cbar_cost / double(n_pairs);
}

/*rank of task*/
double ranku(int id_taski, Data* data, vector<double> &ranku_aux, double lambda) {
    auto f_suc = data->succ.find(id_taski);

    auto rank = [&](int id_taskj) {
        return cbar(id_taski, id_taskj, data, lambda) + ranku(id_taskj, data, ranku_aux, lambda);
    };

    if (f_suc != data->succ.end() && f_suc->second.size() != 0) {
        auto max_value = 0.0;
        for_each(f_suc->second.begin(), f_suc->second.end(), [&](int id_taskj) {
            double val = 0.0;
            if (ranku_aux[id_taskj] == -1)
                val = ranku_aux[id_taskj] = rank(id_taskj);
            else
                val = ranku_aux[id_taskj];

            max_value = std::max(max_value, val);

        });
        // Check if id_taski is root (this ensures that the root task has the greatest rank)
        id_taski == data->id_root ? max_value *= 2 : max_value;
        return wbar(id_taski, data) + max_value;
    } else {
        return wbar(id_taski, data);
    }
}

double find_first_gap(vector<Event> vm_orders, double desired_start_time, double duration) {
    /*Find the first gap in an agent's list of jobs
    The gap must be after `desired_start_time` and of length at least
    duration.
    */

    // No task: can fit it in whenever the job is ready to run
    if (vm_orders.size() == 0)
        return desired_start_time;

    /* Try to fit it in between each pair of Events, but first prepend a
    dummy Event which ends at time 0 to check for gaps before any real
    Event starts.*/
    vector<Event> aux(1);
    auto a = boost::join(aux, vm_orders);

    for (unsigned i = 0; i < a.size() - 1; i++) {
        auto earlist_start = std::max(desired_start_time, a[i].end);
        if (a[i + 1].start - earlist_start > duration)
            return earlist_start;
    }

    // No gaps found: put it at the end, or whenever the task is ready
    return std::max(vm_orders.back().end, desired_start_time);
}

/* Earliest time that task can be executed on vm*/
double start_time(int id_task, int id_vm, vector<int> taskOn, event_map orders, vector<double> end_time, Data* data,
                  double lambda) {
    auto duration = compcost(id_task, id_vm, data);
    auto comm_ready = 0.0;
    auto max_value = 0.0;
    if (id_task != data->id_root && data->prec.find(id_task)->second.size() > 0) {
        //for each prec of task
        for_each(data->prec.find(id_task)->second.begin(), data->prec.find(id_task)->second.end(), [&](const int &p) {
            //comm_ready = std::max(end_time(p, orders.find(taskOn[p])->second) + commcost(p, id_task, taskOn[p], id_vm), comm_ready);
            max_value = std::max(end_time[p], max_value);
            comm_ready += commcost_dynamic(p, id_task, taskOn[p], id_vm, data, lambda);
        });
    }

    auto f = orders.find(id_vm);
    auto queue_value = 0.0;
    if (f != orders.end() && !f->second.empty())
        queue_value = (f->second.back().end);

    max_value = std::max(max_value, queue_value);

    comm_ready += commcost_static(id_task, id_vm, data, lambda);
    comm_ready = comm_ready + max_value;
    return find_first_gap(orders.find(id_vm)->second, comm_ready, duration);
}

/*
 * Allocate task to the vm with earliest finish time
 */
void
allocate(int id_task, vector<int> &taskOn, vector<int> vm_keys, event_map &orders, vector<double> &end_time, Data* data,
         double lambda) {
    auto st = [&](int id_vm) {
        return start_time(id_task, id_vm, taskOn, orders, end_time, data, lambda);
    };
    auto ft = [&](int id_vm) {
        return st(id_vm) + compcost(id_task, id_vm, data);
    };
    //sort vms based on task finish time

    sort(vm_keys.begin(), vm_keys.end(),
         [&](const int &vma, const int &vmb) {
             return ft(vma) < ft(vmb);
         });
    auto vm_id = vm_keys.front();
    auto start = st(vm_id);
    auto end = ft(vm_id);


    Event event;
    event.id = id_task;
    event.start = start;
    event.end = end;

    end_time[id_task] = end;

    auto f = orders.find(vm_id);
    f->second.push_back(event);

    sort(f->second.begin(), f->second.end(), [&](const Event &eventa, const Event &eventb) {
        return eventa.start < eventb.start;
    });
    taskOn[id_task] = vm_id;
}

// Get the next task based on the start time and remove the task
// if there is no task, return -1
int get_next_task(event_map &orders) {

    auto min_start_time = numeric_limits<double>::max();
    int task_id = -1, vm_id;
    for (auto info : orders) {
        if (!info.second.empty()) {
            if (info.second.begin()->start < min_start_time) {
                min_start_time = info.second.begin()->start;
                task_id = info.second.begin()->id;
                vm_id = info.first;
            }
        }
    }
    if (task_id != -1)
        orders.find(vm_id)->second.erase(orders.find(vm_id)->second.begin());

    return task_id;
}


#endif //HEA_HEFT_H
