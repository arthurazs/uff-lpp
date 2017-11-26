/*
 * Data.h
 *
 *	Read the input files (vcl and dag) to building the basic data structures
 *
 *  Created on: Jun 2, 2016
 *      Author: Luan Teylo
 */

#ifndef DATA_H_
#define DATA_H_

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <unordered_map>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/bimap.hpp>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <map>

using namespace std;


/*File object*/
class File{
public:
	string name;
	int id;
	double size;
	bool is_static;
	int static_vm;
    vector<int> all_tasks;

	File(string name, int id, double size, bool is_static = false, int static_vm = -1)
	: name(name), id(id), size(size), is_static(is_static), static_vm(static_vm){}
};

/*Task Object*/
class Task{
public:
	string name;
	int id;
	string tag;
	double base_time;
	vector<int> input, output;
	Task (string name, int id, string tag, double base_time, vector<int> input, vector<int> output):
		name(name), id(id), tag(tag), base_time(base_time), input(input), output(output){}
};

/*Vm object*/
class VMachine{
public:
	string name;
	int id;
	double slowdown, storage, cost, bandwidth;
	int type_id;
	VMachine (string name, int id, double slowdown, double storage, double cost, double bandwidth, int type_id):
		name(name), id(id), slowdown(slowdown), storage(storage), cost(cost), bandwidth(bandwidth), type_id(type_id){}
	//VMachine(VMachine & other):
	//	Obj(other.name, other.id), slowdown(other.slowdown), storage(other.storage), cost(other.cost), bandwidth(other.bandwidth){}
};


class Data {

public:
	unordered_map <int, vector<int> > succ, prec; //Workflow task Graphs
	int task_size, sfile_size, dfile_size, file_size, size, vm_size, id_sink, id_root;
	double period_hr;

	unordered_map<int, File> file_map;
	unordered_map<int, Task> task_map;
	unordered_map<int, VMachine> vm_map;

	unordered_map <string, int> key_map;

	vector<int> static_vms;
	vector<int> height;
    vector<double> storage_vet; //storage of vm

	Data(string* dag_file, string* vcl_file){
		loadFiles(dag_file, vcl_file);
		height.resize(task_size, -1);
		computeHeight(this->id_root, 0);
	}
	virtual ~Data(){}

private:

	/* Read Input Files */
	void loadFiles(string* workflow, string* cluster){
		double total_storage = 0.0;
		double total_file = 0.0;

		//Reading file
		ifstream in_file(*workflow);
		string line;

		// Get number of tasks and number of files
		getline(in_file, line);
		vector<string> tokens;
		boost::split(tokens, line, boost::is_any_of(" "));

		sfile_size = stoi(tokens[0]);
		dfile_size = stoi(tokens[1]);
		task_size = stoi(tokens[2]) + 2;
		file_size = sfile_size + dfile_size;
		size = task_size + dfile_size;

		getline(in_file, line); //reading blank line

		//start initial integer_id of elements
		int id_task = 1;
		int id_dfile = task_size;
		int id_sfile = task_size + dfile_size;

		//Reading files
		for(int i = 0; i < file_size; i++){
			getline(in_file, line);
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" "));
			auto file_name = strs[0];
			auto file_size = stod(strs[1]);

			int id = 0;
			bool is_static = false;
			auto place = -1;

			if(i < sfile_size){
				place = stoi(strs[3]);
				id = id_sfile;
				id_sfile += 1;
				is_static = true;
			}else{
				id = id_dfile;
				id_dfile += 1;
			}

			auto afile = File(file_name, id, file_size, is_static, place);
			total_file += file_size;
			key_map.insert(make_pair(file_name, id));
			file_map.insert(make_pair(id, afile));
		}
		getline(in_file, line); //reading blank line

		//read tasks
		//cout << "Reading tasks" << endl;
		for(int i = 0; i < task_size-2; i++){
			getline(in_file, line);
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" "));
			// get task info
			auto tag = strs[0];
			string task_name = strs[1];
			auto base_time = stod(strs[2]);
			auto in_size = stoi(strs[3]);
			auto out_size = stoi(strs[4]);
			int current_task = id_task;
			id_task += 1;
			//cout << taskId << " " << taskName << " " << taskTime << " " << inSize << " " << outSize << endl;
			//Add task in bi-map idToFile
			vector<int> input;
			vector<int> output;

			//reading input files
			for(int j = 0; j < in_size; j++){
				getline(in_file, line);
				auto fileKey = key_map.find(line)->second;
				input.push_back(fileKey);
			}

			//reading output files
			for(int j = 0; j < out_size; j++){
				getline(in_file, line);
				auto fileKey = key_map.find(line)->second;
				//update file
				output.push_back(fileKey);
			}


			Task atask(task_name, current_task, tag, base_time, input, output);

			key_map.insert(make_pair(tag, current_task));
			task_map.insert(make_pair(current_task, atask));
		}
		getline(in_file, line); //reading blank line

		//Update Root and Sink tasks
		id_root = 0;
		id_sink = id_task;
		id_task += 1;

		key_map.insert(make_pair("root", id_root));
		key_map.insert(make_pair("sink", id_sink));


		Task rootTask("root", id_root, "ROOT", 0.0, vector<int>(0), vector<int>(0));
		Task sinkTask("sink", id_sink, "SINK", 0.0, vector<int>(0), vector<int>(0));

		task_map.insert(make_pair(id_root, rootTask));
		task_map.insert(make_pair(id_sink, sinkTask));

		auto f_root = succ.insert(make_pair(id_root, vector<int>()));
		succ.insert(make_pair(id_sink, vector<int>()));

		vector<int> aux(task_size,-1);
		aux[id_root] = 0;
		aux[id_sink] = 0;

		//reading succ graph
		for(int i = 0; i < task_size-2; i++){
			getline(in_file, line);//Reading parent task
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" "));
			auto task_tag = strs[0];
			auto csize = stoi(strs[1]);

			auto task_id = key_map.find(task_tag)->second;

			vector<int> children;
			//Reading children task
			for(int j = 0; j < csize; j++){
				getline(in_file, line);
				auto child_id = key_map.find(line)->second;
				children.push_back(child_id);
				aux[child_id] = 0;
			}

			//sink task
			if(csize == 0)
				children.push_back(id_sink);
			succ.insert(make_pair(task_id,children));
		}

		//add synthetic root task
		for(int i = 0; i < task_size; i++)
			//add root
			if(aux[i] == -1)
				f_root.first->second.push_back(i);

		prec = reverse_map(succ);
		in_file.close();

		ifstream in_cfile(*cluster);



		// Reading Cluster's info
		getline(in_cfile, line);//ignore first line
		getline(in_cfile, line);

		vector<string> strs1;
		boost::split(strs1, line, boost::is_any_of(" "));

		period_hr = stod(strs1[2]);

		vm_size = stoi(strs1[4]);

		storage_vet.resize(vm_size, 0);

		int vm_id = 0;
		//reading vms
		for(auto i = 0; i < vm_size; i++){
			getline(in_cfile, line);
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" "));


			int type_id = stoi(strs[0]);
			string vm_name = strs[1];
			double slowdown = stod(strs[2]);
			double storage = stod(strs[3]) * 1024; // GB to MB
			double bandwidth = stod(strs[4]);
			double cost = stod(strs[5]);

			VMachine avm(vm_name, vm_id, slowdown, storage, cost, bandwidth, type_id);
			vm_map.insert(make_pair(vm_id, avm));
			storage_vet[vm_id] = storage;
			vm_id += 1;
			total_storage += storage;
		}

		//Check if storage is enough
		for(auto it : file_map){
			if(it.second.is_static){
				storage_vet[it.second.static_vm] -= it.second.size;
				if(storage_vet[it.second.static_vm] < 0){
					cerr << "Static file is bigger than the vm capacity" << endl;
					throw;
				}
			}
		}

		if(total_storage < total_file){
			cerr << "Storage is not enough" << endl;
			throw;
		}

	}

	void computeHeight(int node, int n){
		if(height[node] < n){
			height[node] = n;
			auto vet = succ.find(node)->second;
			for(auto j : vet)
				computeHeight(j, n+1);
		}
	}
	unordered_map<int, vector<int>> reverse_map(unordered_map<int, vector<int>> amap){
		unordered_map<int, vector<int>> r_map;
		for(auto key : amap){
			for(auto val : amap.find(key.first)->second){
				auto f = r_map.insert(make_pair(val, vector<int>(0)));
				f.first->second.push_back(key.first);
			}
		}
		return r_map;
	}
};



#endif /* DATA_H_ */











