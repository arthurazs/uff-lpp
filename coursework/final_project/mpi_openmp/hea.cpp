#include <iostream>
#include <fstream>
#include "HEFT.h"
#include "MinMin.h"
#include <mpi.h>
#include <omp.h>

#define TAGEND 7
#define TAGBEST 8
#define TAGCROMO 9
#define TAGREQ 13
#define BUFFERSIZE 10000
#define NTHREADS 2
#define MAX_TIME 1200 // 20 minutes

void printResult(string arquivo, double fitness, int found, double time_found, int total, double time, int num_chromosomes, int num_elite_set, int num_generations, int foundMyGen, double foundMyTime, int numProcs, int numThreads);
int foundSeqGen = -1; // qual geracao ele encontrou o best
clock_t foundSeqTime; //tempo q encontrou o best
// int iTotal = 0; //quantidade de geracoes executadas
int foundMyGen = 0; //Iteracao que foi encontrado o melhor best
clock_t foundMyTime; //Tempo que foi encontrado o melhor best
double fitnessSeq = 0; //Melhor best do sequencial

int global_gen_count = 0;
clock_t begin_time;

typedef vector<Chromosome> vect_chrom_type;
Data* data;

struct Settings_struct {
    int num_chromosomes; // Number of chromosomes
    int num_generations; // Number of generations
    int num_elite_set; // Max size of Elite-set

    float mutation_probability; // Probability of mutation (for each gene)
    float elitism_rate; // Rate of generated solutions
    int howMany_elistism;
    float alpha; // percentage of chromosomes send to local search procedures
    float localSearch_probability; // Probability of local search

    float time_limit; // Run time limite in seconds
    int print_gen; // What generation will be printed
    bool verbose, start_heuristic; // verbose = print type, start_heuristic = Initial Population type

    long seed; // Random Seed

    double delta = 0.0; // acceptance criteria for Elite-Set (based on distance)
    double lambda; // read and write constant

    Settings_struct() {
        // Default Settings
        // num_chromosomes = NUM_CHRO;
        // num_generations = 100;
        // num_elite_set = num_chromosomes / 2;

        mutation_probability = 0.10;
        elitism_rate = 0.10;
        alpha = 0.30; // Local Search Rate
        localSearch_probability = 0.50;
        time_limit = 7200; // time limit 30 minutes

        print_gen = 10;

        verbose = false;
        start_heuristic = true;

        //howMany_elistism = (int) ceil(num_chromosomes * elitism_rate);

        lambda = 0.000;
    }

};


Settings_struct *setting;

/*  Call HEFT */
Chromosome HEFT(Data* data) {
    //orders
    event_map orders;

    //building and ordering the seqOfTasks
    vector<int> seqOftasks;
    boost::copy(data->task_map | boost::adaptors::map_keys, std::back_inserter(seqOftasks));

    vector<double> ranku_vet(seqOftasks.size(), 0.0);
    vector<double> ranku_aux(seqOftasks.size(), -1);

    for_each(seqOftasks.begin(), seqOftasks.end(), [&](const int &idA) {
        ranku_vet[idA] = ranku(idA, data, ranku_aux, setting->lambda);
    });

    sort(seqOftasks.begin(), seqOftasks.end(), [&](const int &idA, const int &idB) {
        return ranku_vet[idA] < ranku_vet[idB];
    });

    //get all vm keys
    vector<int> vm_keys;
    boost::copy(data->vm_map | boost::adaptors::map_keys, std::back_inserter(vm_keys));

    //build orders struct (event_map)
    for (auto vm_key : vm_keys)
        orders.insert(make_pair(vm_key, vector<Event>()));


    vector<int> taskOn(data->task_size, -1);
    vector<double> end_time(data->task_size, 0);
    for (auto id_task = seqOftasks.rbegin(); id_task != seqOftasks.rend(); id_task++) { // reverse vector
        allocate(*id_task, taskOn, vm_keys, orders, end_time, data, setting->lambda);
    }


    // == build chromosome == //

    Chromosome heft_chr(data, setting->lambda);

    // build allocation
    for (auto info : orders) {
        auto id_vm = info.first;
        for (auto event : info.second) {
            auto task = data->task_map.find(event.id)->second;
            heft_chr.allocation[task.id] = id_vm;
            // update output files;
            for (auto out : task.output)
                heft_chr.allocation[out] = id_vm;
        }
    }


    // build ordering
    heft_chr.ordering.clear();
    // add root
    heft_chr.ordering.push_back(data->id_root);
    int task_id = -1;
    do {
        task_id = get_next_task(orders);
        if (task_id != -1 && task_id != data->id_root && task_id != data->id_sink)
            heft_chr.ordering.push_back(task_id);
    } while (task_id != -1);
    // add sink
    heft_chr.ordering.push_back(data->id_sink);

    heft_chr.computeFitness(true, true);

    return heft_chr;
}

/* Call MinMin */
Chromosome minMinHeuristic(Data* data) {
    list<int> task_list;
    // start task list
    for (auto info : data->task_map)
        task_list.push_back(info.second.id);
    task_list.sort([&](const int &a, const int &b) {
        return data->height[a] < data->height[b]; });

    list<int> avail_tasks;

    vector<double> ft_vector(data->size, 0);
    vector<double> queue(data->vm_size, 0);
    vector<int> file_place(data->size, 0);
    list<int> task_ordering(0);


    //the task_list is sorted by the height(t). While task_list is not empty do
    while (!task_list.empty()) {
        auto task = task_list.front(); //get the first task
        avail_tasks.clear();
        while (!task_list.empty() && data->height[task] == data->height[task_list.front()]) {
            //build list of ready tasks, that is the tasks which the predecessor was finish
            avail_tasks.push_back(task_list.front());
            task_list.pop_front();
        }

        schedule(data, avail_tasks, ft_vector, queue, file_place, task_ordering,
                setting->lambda); //Schedule the ready tasks
    }

    Chromosome minMin_chrom(data, setting->lambda);

    for (int i = 0; i < data->size; i++)
        minMin_chrom.allocation[i] = file_place[i];
    minMin_chrom.ordering.clear();

    minMin_chrom.ordering.insert(minMin_chrom.ordering.end(), task_ordering.begin(), task_ordering.end());
    minMin_chrom.computeFitness(true, true);

    //if(setting->verbose)
    //    cout << "MinMIn fitness: " << minMin_chrom.fitness << endl;

    return minMin_chrom;
}




// ========== Path Relinking ============ //

Chromosome pathRelinking(vect_chrom_type& Elite_set, const Chromosome &dest, Data* data) {

    Chromosome best(dest);

    // For each chromosome on Elite Set, do:
    #pragma omp parallel for
    for (unsigned i = 0; i < Elite_set.size(); i++) {
        auto src = Elite_set[i];

        // Copy ordering from dest chromosome
        src.ordering.clear();
        src.ordering.insert(src.ordering.end(), dest.ordering.begin(), dest.ordering.end());

        for (int el = 0; el < data->size; el++) {
            if (src.allocation[el] != dest.allocation[el]) {
                src.allocation[el] = dest.allocation[el];
                src.computeFitness(true, true);
                #pragma omp critical
                {
                    if (best.fitness > src.fitness)
                        best = src;
                }
            }
        }
    }


    return best;
}

// Get the best chromosome

inline int getBest(vect_chrom_type& Population) {
    Chromosome* best = &(Population[0]);
    auto pos = 0;
    for (int i = 0; i < setting->num_chromosomes; i++)
        if (best->fitness > Population[i].fitness) {
            best = &Population[i];
            pos = i;
        }
    return pos;
}

// Tournament Selection

inline int tournamentSelection(vect_chrom_type& Population) {
    //we pick to chromosomes at random
    int a = random() % Population.size();
    int b = random() % Population.size();

    //make sure they're not the same chromosome!
    while (b == a)
        b = random() % Population.size();

    //now select the better of the two as our parent
    return Population[a].fitness < Population[b].fitness ? a : b;
}


// =========== Local search functions  ========= //

// N1 - Swap-vm

inline Chromosome localSearchN1(const Data* data, Chromosome* ch) {
    Chromosome old_ch(ch);
    for (int i = 0; i < data->size; i++) {
        for (int j = i + 1; j < data->size; j++) {
            if (ch->allocation[i] != ch->allocation[j]) {
                //do the swap
                iter_swap(ch->allocation.begin() + i, ch->allocation.begin() + j);
                ch->computeFitness();
                if (ch->fitness < old_ch.fitness) {
                    return ch;
                }
                //return elements
                iter_swap(ch->allocation.begin() + i, ch->allocation.begin() + j);
            }
        }
    }
    return old_ch;
}

// N2 - Swap position

inline Chromosome localSearchN2(const Data* data, Chromosome* ch) {
    Chromosome old_ch(ch);
    // for each task, do
    for (int i = 0; i < data->task_size; i++) {
        auto task_i = ch->ordering[i];
        for (int j = i + 1; j < data->task_size; j++) {
            auto task_j = ch->ordering[j];
            if (ch->height_soft[task_i] == ch->height_soft[task_j]) {
                //do the swap
                iter_swap(ch->ordering.begin() + i, ch->ordering.begin() + j);
                ch->computeFitness(false, true);
                if (ch->fitness < old_ch.fitness) {
                    return ch;
                }
                //return elements
                iter_swap(ch->ordering.begin() + i, ch->ordering.begin() + j);
            } else
                break;
        }
    }
    return old_ch;
}

// N3 = Move-1 Element

inline Chromosome localSearchN3(const Data* data, Chromosome* ch) {
    Chromosome old_ch(ch);

    for (int i = 0; i < data->size; ++i) {
        int old_vm = ch->allocation[i];
        for (int j = 0; j < data->vm_size; j++) {
            if (old_vm != j) {
                ch->allocation[i] = j;
                ch->computeFitness();
                if (ch->fitness < old_ch.fitness) {
                    return ch;
                }
            }
        }
        ch->allocation[i] = old_vm;
    }
    return old_ch;
}



// ========== Main Functions ========== //

inline void doNextPopulation(vect_chrom_type &Population) {
    //int how_many =  (int) ceil(setting->num_chromosomes * setting->elitism_rate);

    vector<Chromosome> children_pool;

    // === do offsprings === //
    for (int i = 0; i < ceil(setting->num_chromosomes / 2.0); i++) {
        // select our two parents with tournament Selection
        int posA, posB;
        posA = tournamentSelection(Population);
        do {
            posB = tournamentSelection(Population);
        } while (posA == posB);

        // get the parents
        auto parentA = Population[posA];
        auto parentB = Population[posB];

        // cross their genes
        auto child = parentA.crossover(&parentB);
        // mutate the child
        child.mutate(setting->mutation_probability);
        // recompute fitness
        child.computeFitness();
        // Add solution on children_pool
        children_pool.push_back(child);
    }

    // === update population === //

    // add all solutions to the children_pool
    children_pool.insert(children_pool.end(), Population.begin(), Population.end());

    // Delete old population
    Population.clear();

    // Elitisme operator - the best is always on the population
    //auto posBest = getBest(children_pool);
    //Population.push_back(children_pool[posBest]);
    sort(children_pool.begin(), children_pool.end(), [&](const Chromosome &chr1, const Chromosome & chr2) {
        return chr1.fitness < chr2.fitness;
    });

    for (int i = 0; i < setting->howMany_elistism; i++) {
        Population.push_back(children_pool[0]);
        children_pool.erase(children_pool.begin());
    }

    // Selected the solutions to build the new population
    while (Population.size() < static_cast<unsigned int> (setting->num_chromosomes)) {
        auto pos = tournamentSelection(children_pool);
        Population.push_back(Chromosome(children_pool[pos]));
        children_pool.erase(children_pool.begin() + pos);
    }
    random_shuffle(Population.begin(), Population.end());
}

// Call all Local Search Functions

inline void localSearch(vect_chrom_type &Population, Data* data) {

    int how_many = setting->alpha * setting->num_chromosomes;

    for (int j = 0; j < how_many; j++) {
        auto ch_pos = tournamentSelection(Population);
        Population[ch_pos] = localSearchN1(data, &Population[ch_pos]);
        Population[ch_pos] = localSearchN2(data, &Population[ch_pos]);
        Population[ch_pos] = localSearchN3(data, &Population[ch_pos]);
    }
}

Chromosome run(string* name_workflow, string* name_cluster, int my_rank, int num_procs) {

    // Load input Files and the data structures used by the algorithms
    data = new Data(name_workflow, name_cluster);

    vector<Chromosome> Population;
    int aux_size = (ceil(setting->num_chromosomes * 0.9) * 2) + ceil(setting->num_chromosomes * 0.1) + 2;
    Population.resize(aux_size);

    vector<Chromosome> Elite_set;
    int flag;

    // Set Delta
    setting->delta = data->size / 4.0;

    // check distance (inner Function)
    auto check_distance = [&](Chromosome* chr, const vector<Chromosome>* Set) {
        for (auto set_ch : *Set) {
            if (chr->getDistance(set_ch) < setting->delta) {
                return false;
            }

        }
        return true;
    };

    // == Start initial population == //
    Chromosome minminChr(minMinHeuristic(data));
    Chromosome heftChr(HEFT(data));


    #pragma omp parallel sections
    {
        #pragma omp section
        {
            Population[0] = minminChr;
            Population[1] = heftChr;
        }
        #pragma omp section
        {
            int aux_first = 2;
            int aux_last = (ceil(setting->num_chromosomes * 0.9) * 2) + aux_first;

            double mut = 0.05;

            // 90% using the mutate procedure with variation
            for (int i = aux_first; i < aux_last; i+=2) {
                Chromosome chr1(minminChr);
                chr1.mutate(mut);
                chr1.computeFitness();

                Chromosome chr2(heftChr);
                chr2.mutate(mut);
                chr2.computeFitness();

                Population[i] = chr1;
                Population[i+1] = chr2;

                mut = mut > 1.0 ? 0.05 : mut + 0.05;
            }
        }
        #pragma omp section
        {
            int aux_first = (ceil(setting->num_chromosomes * 0.9) * 2) + 2;
            int aux_last = ceil(setting->num_chromosomes * 0.10) + aux_first;

            // 10% random solutions
            for (int i = aux_first; i < aux_last; i++) {
                Population[i] = Chromosome(data, setting->lambda);
            }
        }
    }

    // Get best solution from initial population
    Chromosome best(Population[getBest(Population)]);

    // Do generation
    int count = 0;
    // start stop clock


    int myEnd = 0;
    int netEnd = 0;
    double local_best;
    int x = 10; //porcetagem de melhora do fitness
    int interval = 25; //intervalo para atualização periodica

    MPI_Status status;


    //sai do while somente quando todas as threads de todos processos
    //atingirem o numero de geracoes esperado
    while (netEnd < num_procs) {

        local_best = best.fitness;

        // Do local Search ?

        float doit = (float) random() / (float) RAND_MAX;

        if (doit <= (setting->localSearch_probability))
            localSearch(Population, data);

        // Update best
        auto pos = getBest(Population);

        if (best.fitness > Population[pos].fitness) {

            best = Population[pos];

            // Apply path Relinking
            if (!Elite_set.empty())
                best = pathRelinking(Elite_set, best, data);

            // Update Elite-set
            if (check_distance(&best, &Elite_set))
                Elite_set.push_back(best); // Push all best' solutions on Elite-set

            // check elite set size
            if (Elite_set.size() > static_cast<unsigned int> (setting->num_elite_set))
                Elite_set.erase(Elite_set.begin());


            // Apply Local Search

            best = localSearchN1(data, &best);
            best = localSearchN2(data, &best);
            best = localSearchN3(data, &best);

            Population[pos] = best;

            count = 0;

            // my best
            foundMyGen = global_gen_count;
            foundMyTime = clock();
        }

        if (best.fitness <= fitnessSeq && foundSeqGen == -1) {
            foundSeqGen = global_gen_count;
            foundSeqTime = clock();
        }

        //verifica se terminou e indica o fim do processo
        double time_running = double(clock() - begin_time) / CLOCKS_PER_SEC;
        if (count == setting->num_generations && myEnd == 0 || time_running > MAX_TIME && myEnd == 0) {
            myEnd++;
            netEnd++;
            for (int i = 0; i < num_procs; i++) {
                if (i != my_rank) {
                    MPI_Send(&myEnd, 1, MPI_INT, i, TAGEND, MPI_COMM_WORLD);

                }
            }
            if (setting->verbose) {
                printf(
                        "|P%d| Generations completed\n",
                        my_rank);
            }
        }

        //troca de mensagens entre processos
        int position;
        char buffer[BUFFERSIZE];

        //verifica se tem mensagem de fim de outros processos e recebe
        do {
            MPI_Iprobe(MPI_ANY_SOURCE, TAGEND, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                int temp;
                MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, TAGEND, MPI_COMM_WORLD, &status);
                netEnd += temp;
            }
        } while (flag);

        //envia o best a cada N iterações ou se a solução melhorar > x%
        if ((global_gen_count % interval == 0) || (((best.fitness / 100) * x) < (local_best - best.fitness))) {
            for (int i = 0; i < num_procs; i++) {
                if (i != my_rank) {
                    MPI_Send(&local_best, 1, MPI_DOUBLE, i, TAGBEST, MPI_COMM_WORLD);
                }
            }
        }

        //verifica se tem mensagem de fitness e recebe
        //se o fitness recebido for melhor, pede o cromossomo
        do {
            MPI_Iprobe(MPI_ANY_SOURCE, TAGBEST, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                double aux_best;
                MPI_Recv(&aux_best, 1, MPI_DOUBLE, MPI_ANY_SOURCE, TAGBEST, MPI_COMM_WORLD, &status);
                if (aux_best < best.fitness) {
                    int aux = 0;
                    MPI_Send(&aux, 1, MPI_INT, status.MPI_SOURCE, TAGREQ, MPI_COMM_WORLD);
                }
            }
        } while (flag);

        //verifica se tem pedido de cromossomo e envia
        do {
            MPI_Iprobe(MPI_ANY_SOURCE, TAGREQ, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                int aux;
                MPI_Recv(&aux, 1, MPI_INT, MPI_ANY_SOURCE, TAGREQ, MPI_COMM_WORLD, &status);

                //packing
                position = 0;

                // double
                MPI_Pack(&best.fitness, 1, MPI_DOUBLE, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);
                MPI_Pack(&best.lambda, 1, MPI_DOUBLE, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);
                MPI_Pack(&best.transfer_size, 1, MPI_DOUBLE, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);

                // int vector
                int alloc_size = best.allocation.size();
                int height_size = best.height_soft.size();
                int order_size = best.ordering.size();
                MPI_Pack(&alloc_size, 1, MPI_INT, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);
                MPI_Pack(best.allocation.data(), alloc_size, MPI_INT, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);
                MPI_Pack(&height_size, 1, MPI_INT, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);
                MPI_Pack(best.height_soft.data(), height_size, MPI_INT, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);
                MPI_Pack(&order_size, 1, MPI_INT, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);
                MPI_Pack(best.ordering.data(), order_size, MPI_INT, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);

                // double vector
                int time_size = best.time_vector.size();
                int start_size = best.start_time_vector.size();
                MPI_Pack(&time_size, 1, MPI_INT, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);
                MPI_Pack(best.time_vector.data(), time_size, MPI_INT, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);
                MPI_Pack(&start_size, 1, MPI_INT, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);
                MPI_Pack(best.start_time_vector.data(), start_size, MPI_INT, buffer, BUFFERSIZE, &position, MPI_COMM_WORLD);

                // send pack
                MPI_Send(buffer, BUFFERSIZE, MPI_PACKED, status.MPI_SOURCE, TAGCROMO, MPI_COMM_WORLD);
            }
        } while (flag);

        //TODO
        do {
            MPI_Iprobe(MPI_ANY_SOURCE, TAGCROMO, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                //packing
                position = 0;
                Chromosome aux;
                MPI_Recv(&buffer, BUFFERSIZE, MPI_PACKED, MPI_ANY_SOURCE, TAGCROMO, MPI_COMM_WORLD, &status);

                // double
                MPI_Unpack(buffer, BUFFERSIZE, &position, &aux.fitness, 1, MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Unpack(buffer, BUFFERSIZE, &position, &aux.lambda, 1, MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Unpack(buffer, BUFFERSIZE, &position, &aux.transfer_size, 1, MPI_DOUBLE, MPI_COMM_WORLD);

                // int vector
                int alloc_size;
                int height_size;
                int order_size;
                MPI_Unpack(buffer, BUFFERSIZE, &position, &alloc_size, 1, MPI_INT, MPI_COMM_WORLD);
                aux.allocation.resize(alloc_size);
                MPI_Unpack(buffer, BUFFERSIZE, &position, &aux.allocation[0], alloc_size, MPI_INT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, BUFFERSIZE, &position, &height_size, 1, MPI_INT, MPI_COMM_WORLD);
                aux.height_soft.resize(height_size);
                MPI_Unpack(buffer, BUFFERSIZE, &position, &aux.height_soft[0], height_size, MPI_INT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, BUFFERSIZE, &position, &order_size, 1, MPI_INT, MPI_COMM_WORLD);
                aux.ordering.resize(order_size);
                MPI_Unpack(buffer, BUFFERSIZE, &position, &aux.ordering[0], order_size, MPI_INT, MPI_COMM_WORLD);

                // double vector
                int time_size;
                int start_size;
                MPI_Unpack(buffer, BUFFERSIZE, &position, &time_size, 1, MPI_INT, MPI_COMM_WORLD);
                aux.time_vector.resize(time_size);
                MPI_Unpack(buffer, BUFFERSIZE, &position, &aux.time_vector[0], time_size, MPI_INT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, BUFFERSIZE, &position, &start_size, 1, MPI_INT, MPI_COMM_WORLD);
                aux.start_time_vector.resize(start_size);
                MPI_Unpack(buffer, BUFFERSIZE, &position, &aux.start_time_vector[0], start_size, MPI_INT, MPI_COMM_WORLD);

                if (best.fitness > aux.fitness) {
                    Population.push_back(aux);
                    auto pos = getBest(Population);
                    best = Population[pos];
                }
            }
        } while (flag);


        if (setting->verbose && (count % setting->print_gen) == 0) {
            cout << "|P" << my_rank << "|" << "Gen: " << count << " Fitness: " << best.fitness / 60.0 << "(s)" << endl;
        }

        count++;
        global_gen_count++;

        doNextPopulation(Population);

    }
    if (setting->verbose) {
        printf(
                "\n|P%d| Stopping\n",
                my_rank);
    }

    // return the global best
    return best;
}




// Read command line parameters (input files)

void setupCmd(int argc, char **argv, string* name_workflow, string* name_cluster) {

    try {

        // Define the command line object.
        CmdLine cmd("Hybrid Evolutionary Algorithm", ' ', "1.0");

        // Define a value argument and add it to the command line.
        ValueArg<string> arg1("w", "workflow", "Name of workflow file", true, "file", "string");
        cmd.add(arg1);
        ValueArg<string> arg2("c", "cluster", "Name of virtual cluster file", true, "file", "string");
        cmd.add(arg2);

        ValueArg<string> arg3("p", "Populacao", "Qtd populacao", true, "file", "string");
        cmd.add(arg3);
        ValueArg<string> arg4("g", "Geracoes", "Qtd de Geracoes", true, "file", "string");
        cmd.add(arg4);
        ValueArg<string> arg5("b", "Best", "Best", true, "file", "string");
        cmd.add(arg5);

        SwitchArg verbose_arg("v", "verbose", "Output info", cmd, false);

        // Parse the args.
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        *name_workflow = arg1.getValue();
        *name_cluster = arg2.getValue();
        setting->verbose = verbose_arg.getValue();

        std::string::size_type sz;
        setting->num_chromosomes = std::stoi(arg3.getValue(), &sz);
        setting->num_elite_set = setting->num_chromosomes / 2;
        setting->num_generations = std::stoi(arg4.getValue(), &sz);
        fitnessSeq = std::stod(arg5.getValue(), &sz);

    } catch (ArgException &e) { // catch any exceptions
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    }
}

int main(int argc, char **argv) {

    int numprocs, rank, master = 0, tag = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Status status;
    MPI_Comm_rank(world, &rank);
    MPI_Comm_size(world, &numprocs);

    clock_t begin = clock();
    begin_time = clock();
    foundSeqTime = begin;
    foundMyTime = begin;

    string name_workflow, name_cluster;

    setting = new Settings_struct();

    setupCmd(argc, argv, &name_workflow, &name_cluster);

    setting->num_chromosomes = setting->num_chromosomes / numprocs;
    setting->num_elite_set = setting->num_chromosomes / 2;
    setting->howMany_elistism = (int) ceil(setting->num_chromosomes * setting->elitism_rate);

    Chromosome best = run(&name_workflow, &name_cluster, rank, numprocs);
    best.computeFitness(true, true);

    int printing = false;

    if (rank != master) {
        MPI_Send(&best.fitness, 1, MPI_DOUBLE, master, tag, world);
    } else {

        // Compare all fitnesses
        double fitnesses[numprocs];
        fitnesses[master] = best.fitness;
        for (int i = 1; i < numprocs; i++) {
            double local;
            MPI_Recv(&local, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tag, world, &status);
            fitnesses[status.MPI_SOURCE] = local;
        }

        double global_best_fitness = fitnesses[master];
        int global_position = master;
        for (int i = 0; i < numprocs; i++) {
            if (fitnesses[i] < global_best_fitness) {
                global_best_fitness = fitnesses[i];
                global_position = i;
            }
        }

        for (int i = 1; i < numprocs; i++) {
            bool printing_local = false;
            if (i == global_position)
                printing_local = true;
            MPI_Send(&printing_local, 1, MPI_INT, i, tag, world);
        }

        if (global_position == master)
            printing = true;
    }

    if (rank != master) {
        MPI_Recv(&printing, 1, MPI_INT, master, tag, world, &status);
    }

    if (printing == true) {
        clock_t end = clock();

        double elapseSecs = double(end - begin) / CLOCKS_PER_SEC;

        if (setting->verbose) {
            cout << "\t **** HEA **** " << endl;
            best.print();
            printf("\n");
        }

        double elapseTimeFoundSecs = double(foundSeqTime - begin) / CLOCKS_PER_SEC;
        double elapseMyTimeFoundSecs = double(foundMyTime - begin) / CLOCKS_PER_SEC;
        printResult(name_workflow, best.fitness, foundSeqGen, elapseTimeFoundSecs, global_gen_count, elapseSecs, setting->num_chromosomes, setting->num_elite_set, setting->num_generations, foundMyGen, elapseMyTimeFoundSecs, numprocs, NTHREADS);

        if (setting->verbose) {
            cout << "Best fitness found by rank: " << rank << endl;
        }
        cout << "Best fitness: " << best.fitness / 60.0 << "(min)" << " Runtime: " << elapseSecs << "(sec)" << endl;
    }

    delete data;
    //delete setting struct
    delete[] setting;

    MPI_Finalize();
    return 0;
}

void printResult(string arquivo, double fitness, int found, double time_found, int total, double time, int num_chromosomes, int num_elite_set, int num_generations, int foundMyGen, double foundMyTime, int numProcs = 1, int numThreads = 1) {

    std::ofstream out;
    out.open("csv/mpi_openmp.csv", std::ios::app);
    // Arquivo;Processos;Threads;Fitness (Sequencial);Geração (Sequencial);Tempo (Sequencial);Fitness (Meu);Geração (Meu);Tempo (Meu);Gerações (Total);Tempo (Total);Cromossomos;Conjunto Elite;Mínimo de Gerações
    out << arquivo << ";" << numProcs << ';' << numThreads << ';' << fitnessSeq << ";" << found << ';' << time_found << ';' << fitness << ';' << foundMyGen << ";" << foundMyTime << ";" << total << ';' << time << ';' << num_chromosomes << ';' << num_elite_set << ';' << num_generations << "\n";
    out.close();

}
