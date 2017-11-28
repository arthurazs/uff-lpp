/*
* Chrom.h
*
*  Created on: 14 de mar de 2016
*      Author: luan
*/

#ifndef CHROM_H_
#define CHROM_H_

#include <iostream>
#include <vector>
#include <list>
#include <limits>
#include <algorithm>
#include <random>
#include "Data.h"
#include <omp.h>

using namespace std;

struct vm_time {
    vector<pair<double, double>> life;
    double last;
    double static_vm = false;
    double storage = 0.0;
    int count_file = 0;
};

std::random_device rd_chr;
std::mt19937 engine_chr(rd_chr());

class Chromosome {

public:
    vector<int> allocation, height_soft;
    vector<double> time_vector, start_time_vector;
    vector<int> ordering;
    double fitness, lambda, transfer_size;
    Data* data;
    unordered_map<int, vector<string>> scheduler;
    unordered_map<int, vector<int>> vm_queue;


    Chromosome(Data* data, double lambda) :
            allocation(data->size, -1), height_soft(data->task_size, -1),
            time_vector(data->task_size, -1), start_time_vector(data->task_size, -1),
            ordering(0), fitness(0), lambda(lambda), transfer_size(0), data(data) {
        computeHeightSoft(data->id_root);
        encode();
        computeFitness();
    }


    Chromosome(const Chromosome* other) :
            allocation(other->allocation), height_soft(other->height_soft),
            time_vector(other->time_vector), start_time_vector(other->start_time_vector),
            ordering(other->ordering), fitness(other->fitness),
            lambda(other->lambda), transfer_size(other->transfer_size), data(other->data),
            scheduler(other->scheduler), vm_queue(other->vm_queue) {}

    Chromosome() {}


    //virtual ~Chromosome();
    virtual ~Chromosome() {
        // TODO Auto-generated destructor stub
    }

    // Compute the fitness of chromosome
    void computeFitness(bool check_storage = true, bool check_sequence = false) {

        fill(time_vector.begin(), time_vector.end(), -1);
        fill(start_time_vector.begin(), start_time_vector.end(), -1);

        vector<double> queue(data->vm_size, 0);

        if (check_storage && !checkFiles()) {
            std::cerr << "check file error" << endl;
            throw;
        }

        scheduler.clear();
        vm_queue.clear();
        // compute makespan
        for (auto id_task : ordering) {//for each task, do
            if (id_task != data->id_root && id_task != data->id_sink) {//if is not root or sink than

                if (check_sequence && !checkTaskSeq(id_task)) {
                    std::cerr << "Encode error - Chrom: Error in the precedence relations." << endl;
                    throw;
                }

                // Load Vm
                auto vm = data->vm_map.find(allocation[id_task])->second;
                auto task = data->task_map.find(id_task)->second;

                // update vm queue
                auto f_queue = vm_queue.insert(make_pair(vm.id, vector<int>()));
                f_queue.first->second.push_back(task.id);

                // update scheduler
                auto f_scheduler = scheduler.insert(make_pair(vm.id, vector<string>()));
                f_scheduler.first->second.push_back(task.tag);

                // Compute Task Times
                auto start_time = ST(&task, &vm, queue);
                auto read_time = RT(&task, &vm);
                auto run_time = ceil(task.base_time * vm.slowdown);//seconds
                auto write_time = WT(&task, &vm);
                auto finish_time = start_time + read_time + run_time + write_time;

                //update structures
                time_vector[id_task] = finish_time;
                start_time_vector[id_task] = start_time;
                queue[vm.id] = finish_time;

            } else {// root and sink tasks
                if (id_task == data->id_root)
                    time_vector[id_task] = 0;
                else {//sink task
                    double max_value = 0.0;
                    for (auto tk : data->prec.find(id_task)->second)
                        max_value = std::max(max_value, time_vector[tk]);
                    time_vector[id_task] = max_value;
                }
            }
        }

        fitness = time_vector[data->id_sink];
    }

    /*crossover*/
    Chromosome crossover(Chromosome* partner) {
        Chromosome chr(partner);
        uniform_int_distribution<> dis_ordering(0, data->task_size - 1);
        uniform_int_distribution<> dis_allocation(0, data->size - 1);

        int point_ordering = dis_ordering(engine_chr);//crossover point to ordering list
        int point_allocation = dis_allocation(engine_chr);//crossover point to allocation list
        //allocation crossover (Single point crossover)
        for (int i = 0; i < point_allocation; i++) {
            chr.allocation[i] = allocation[i];
        }
        //ordering crossover
        vector<bool> aux(data->task_size, false);
        chr.ordering.clear();
        //ordering crossover first part self -> chr
        for (auto i = 0; i < point_ordering; i++) {
            chr.ordering.push_back(ordering[i]);
            aux[ordering[i]] = true;
        }

        //Ordering crossover second part partner -> chr
        for (auto i = 0; i < data->task_size; i++) {
            if (!aux[partner->ordering[i]])
                chr.ordering.push_back(partner->ordering[i]);
        }
        return chr;
    }

    //Mutation on allocation chromosome
    void mutate(double prob) {
        uniform_int_distribution<> idis(0, data->vm_size - 1);
        for (int i = 0; i < data->size; i++) {
            if (((float) random() / (float) RAND_MAX) <= prob) {
                allocation[i] = idis(engine_chr);
            }
        }
    }


    void print() {

        cout << "#!# " << fitness << endl;
        cout << "Tasks: " << endl;
        for (auto info : data->vm_map) {
            auto vm = info.second;
            cout << vm.id << ": ";
            auto f = scheduler.find(vm.id);
            if (f != scheduler.end()) {
                for (auto task_tag: f->second)
                    cout << task_tag << " ";
            }
            cout << endl;
        }

        cout << endl;

        /*for(auto info : data->vm_map){
            auto vm = info.second;
            cout << "[" << vm.id << "]" << " <" << vm.name << "> : ";
            auto f = scheduler.find(vm.id);
            if(f != scheduler.end()) {
                for (auto task_tag : f->second)
                    cout << task_tag << " ";
            }
            cout << endl;
        }*/
        cout << "Files: " << endl;
        for (auto info: data->vm_map) {
            auto vm = info.second;
            cout << vm.id << ": ";
            for (auto info : data->file_map) {
                auto file = info.second;
                int vm_id = file.is_static ? file.static_vm : allocation[file.id];
                if (vm_id == vm.id)
                    cout << file.name << " ";
            }
            cout << endl;
        }


        /*for(auto info : data->file_map){
            auto file = info.second;
            int vm_id = file.is_static ? file.static_vm : allocation[file.id];
            cout << "["  << vm_id << ", " << file.name << "]" << " ";
        }*/
        //cout << endl;

        /*cout << "Task Sequence: " << endl;
        for(auto task_id : ordering)
            if(task_id != data->id_root && task_id && data->id_sink)
                cout << data->task_map.find(task_id)->second.name <<  ", ";
        cout << endl;*/

    }

    // Compute distance between two solutions
    int getDistance(const Chromosome &chr) {
        int distance = 0;
        #pragma omp parallel sections reduction(+:distance)
        {
            #pragma omp section
            {
                // compute the distance based on position
                for (int i = 0; i < data->size; i++) {
                    if (chr.allocation[i] != allocation[i]) {
                        distance += 1;
                    }
                }
            }
            #pragma omp section
            {
                // compute the distance based on swaps required
                vector<int> aux_ordering(ordering);
                for (int i = 0; i < data->task_size; i++) {
                    if (chr.ordering[i] != aux_ordering[i]) {
                        distance += 1;
                        for (int j = i + 1; j < data->task_size; j++) {
                            if (chr.ordering[i] == aux_ordering[j]) {
                                iter_swap(aux_ordering.begin() + i, aux_ordering.begin() + j);
                            }
                        }
                    }
                }
            }
        }
        return distance;
    }


private:


    int computeHeightSoft(int node) {
        int min = numeric_limits<int>::max();
        if (height_soft[node] != -1)
            return height_soft[node];
        if (node != data->id_sink) {
            for (auto i : data->succ.find(node)->second) {
                int value = computeHeightSoft(i);
                min = std::min(value, min);
            }
        } else {
            height_soft[node] = data->height[node];
            return height_soft[node];
        }
        uniform_int_distribution<> dis(data->height[node], min - 1);
        height_soft[node] = dis(engine_chr);
        return height_soft[node];
    }

    // Random encode (new chromosome)
    void encode() {
        vector<int> seq_list(boost::counting_iterator<int>(0u), boost::counting_iterator<int>(height_soft.size()));
        sort(begin(seq_list), end(seq_list), [&](const int &i, const int &j) {
            return height_soft[i] < height_soft[j];
        });//sort a list based on height soft

        for (auto task : seq_list)//encode ordering Chromosome
            ordering.push_back(task);

        //Encode allocation chromosome
        uniform_int_distribution<> dis(0, data->vm_size - 1);
        for (int i = 0; i < data->size; i++)
            allocation[i] = dis(engine_chr);

    }

    /* Checks the sequence of tasks is valid */
    inline bool checkTaskSeq(int task) {
        for (auto tk : data->prec.find(task)->second)
            if (time_vector[tk] == -1)
                return false;
        return true;
    }

    // Check and organize the file based on the storage capacity
    inline bool checkFiles() {
        bool flag = true;
        int count = 0;

        vector<double> aux_storage(data->storage_vet);
        vector<int> aux(data->vm_size);
        iota(aux.begin(), aux.end(), 0); // 0,1,2,3,4 ... n


        unordered_map<int, vector<int> > map_file;

        int id_vm;
        // build file map and compute the storage
        for (auto it : data->file_map) {
            if (!it.second.is_static) {
                id_vm = allocation[it.second.id];
                auto f = map_file.insert(make_pair(id_vm, vector<int>()));
                f.first->second.push_back(it.second.id);
                auto file = data->file_map.find(it.second.id)->second;
                aux_storage[id_vm] -= file.size;
            }
        }

        do {
            //sort machines based on the storage capacity
            sort(aux.begin(), aux.end(), [&](const int &a, const int &b) {
                return aux_storage[a] < aux_storage[b];
            });

            if (aux_storage[aux[0]] < 0) {//if storage is full, start the file heuristic
                cout << "Starting file heuristic ..." << endl;
                int old_vm = aux[0]; //critical machine
                int new_vm = aux[aux.size() - 1];//best machine

                auto vet_file = map_file.find(old_vm)->second;

                double min = numeric_limits<double>::max();
                int min_file = -1;

                //search min file (based on the size of file)
                for_each(vet_file.begin(), vet_file.end(), [&](int i) {
                    cout << i << endl;
                    auto file = data->file_map.find(i)->second;
                    cout << file.name << endl;
                    if (file.size < min) {
                        min = file.size;
                        min_file = file.id;
                    }
                });

                auto file_min = data->file_map.find(min_file)->second;

                cout << file_min.name << endl;
                //minFile will be move to machine with more empty space
                allocation[file_min.id] = new_vm;
                //update aux Storage
                aux_storage[old_vm] += file_min.size;
                aux_storage[new_vm] -= file_min.size;
                //Update mapFile structure
                map_file[old_vm].erase(remove(map_file[old_vm].begin(), map_file[old_vm].end(), min_file),
                                       map_file[old_vm].end());
                map_file[new_vm].push_back(min_file);
            } else flag = false;

            count++;

        } while (flag && count < data->file_size);

        return !flag;

    }

    // Start Time
    inline double ST(Task* task, VMachine* vm, const vector<double> &queue) {
        // compute wait time
        double max_pred_time = 0.0;

        for (auto tk : data->prec.find(task->id)->second)
            max_pred_time = std::max(max_pred_time, time_vector[tk]);
        return std::max(max_pred_time, queue[vm->id]);
    }

    // Read time
    inline double RT(Task* task, VMachine* vm) {
        //compute read time
        double read_time = 0;
        int id_vm_file;


        for (auto id_file : task->input) {
            auto file = data->file_map.find(id_file)->second;

            if (!file.is_static) {
                id_vm_file = allocation[file.id];
                // update vm queue
                auto f_queue = vm_queue.insert(make_pair(id_vm_file, vector<int>()));
                f_queue.first->second.push_back(file.id);
            } else
                id_vm_file = file.static_vm;


            auto vm_file = data->vm_map.find(id_vm_file)->second;

            read_time += ceil(TT(&file, vm, &vm_file) + (file.size * lambda));
        }

        return read_time;

    }  //Write time

    // Write time
    inline double WT(Task* task, VMachine* vm) {
        //compute the write time
        double write_time = 0;
        for (auto id_file :task->output) {
            auto file = data->file_map.find(id_file)->second;
            auto vm_file = data->vm_map.find(allocation[file.id])->second;

            // update vm queue
            auto f_queue = vm_queue.insert(make_pair(vm_file.id, vector<int>()));
            f_queue.first->second.push_back(id_file);

            write_time += ceil(TT(&file, vm, &vm_file) + (file.size * (lambda * 2)));
        }

        return write_time;
    }

    // Transfer Time
    inline double TT(File* file, VMachine* vm1, VMachine* vm2) {
        if (vm1->id != vm2->id) {// if vm1 != vm2
            //get the smallest link
            auto link = std::min(vm1->bandwidth, vm2->bandwidth);
            return file->size / link;
        }
        return 0;//otherwise
    }


};

#endif /* CHROM_H_ */
