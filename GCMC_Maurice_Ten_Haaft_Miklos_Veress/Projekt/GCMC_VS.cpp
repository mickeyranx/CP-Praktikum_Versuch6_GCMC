

#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <fstream>
#include <list>
#include <bitset>
#include <unordered_map>
#include <omp.h>
#include <string>
using namespace std;


//----------------------------------
//       global params
//----------------------------------

const int M = 64; //lattice length -> size is 64 x 64
static const int particle_width = 1;
static const int particle_length = 8;

//initialize empty OCC lattice
static void init_OCC(bitset<M> (& OCC)[M], int size) {
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++)
        {
            OCC[i][j] = 0;
        }
    }
   
}

//print the OCC lattice to the console in Matrix form
static void print_OCC(bitset<M>(&OCC)[M]) {
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {

            cout << OCC[j][i];
        }

        cout << "\n";
    }
}

//checks the occupancy condition horizontally
static bool checkOccupancyHorizontal(bitset<M> (& OCC)[M], int particle_length ,vector<int> pos) {
    //assign selected lattice point to new variables
    int x = pos[0]; 
    int y = pos[1];
  
    for (int i = 0; i < particle_length; i++)//loop throug the right lattice positions from the chosen point and check wether one of the values is 1
    {
        int current_x = x + i;
        if (current_x > M - 1)  current_x = x + i - M; //jump to opposite lattice point if exceeding lattice limit
        if (OCC[current_x][y] == 1) { 
            return true;
        }
    }
    
    return false; //return false if selected placement ist not occupied
   
    //option 2 should be faster but needs testing
    /*
    for (int i = 0; i < particle_length; i++)
    {
        int current_x = x + i;
        if (current_x > M - 1)  current_x = x + i - M;
        if ((OCC[current_x] >> y & bitset<M>((1 << particle_length) - 1)).count() > 0) {
            return true;
        }
    }
    return false;
    */
}
//checks the occupancy condition vertically
static bool checkOccupancyVertical(bitset<M>(&OCC)[M], int particle_length, vector<int> pos) {
    int x = pos[0];
    int y = pos[1];

    for (int i = 0; i < particle_length; i++)
    {
        int current_y = y + i;
        if (current_y > M - 1) current_y = y + i - M;
        if (OCC[x][current_y] == 1) {
            return true;
        }
    }
    return false;
}


// pos = {x,y}, adds a particle to the OCC lattice with pos being its defining point and return a boolean (successful placement -> true)
static bool add_particle(bitset<M> (&OCC)[M], vector<int> pos, bool horizontal) {
    //assign selected lattice positon to new varables  
    int x = pos[0];
    int y = pos[1];
    if (horizontal && !checkOccupancyHorizontal(OCC, particle_length, pos)) //if horizontal and can be placed
    {
        //->place particle by assigning all its occupying lattice points to 1
        for (int i = 0; i < particle_length; i++)
        {
            int current_x = x + i;
            if (current_x > M - 1) current_x = x + i - M; //jump to opposite lattice point if exceeding lattice limit
            OCC[current_x][y] = 1;
        }
        return true;
    }
    
    else if(!horizontal && !checkOccupancyVertical(OCC, particle_length, pos)){//else if particle is vertical and can be placed
        //->place particle
        for (int i = 0; i < particle_length; i++)
        {
            int current_y = y + i;
            if (current_y > M - 1) current_y = y + i - M; //jump to opposite lattice point if exceeding lattice limit
            OCC[x][current_y] = 1; 
        }
        return true;
    }
    else {
        return false; //return false if particle could not be placed
    }
    
}

//removes a particle from the OCC lattice, prop[i] are the properties of the selected particle {x_pos, y_pos, orientation}
static void remove_particle(bitset<M> (&OCC)[M], vector<int> prop) {
    //assign position properties of selected particle to new variables
    int x = prop[0];
    int y = prop[1];
  
    
    if (prop[2]) //if horizontal 
    {
        //-> remove the particle by changing all its values to 0
        for (int i = 0; i < particle_length; i++)
        {
            int current_x = x + i;
            if (current_x > M - 1) current_x = x + i - M; //jump to opposite lattice point if exceeding lattice limit
            OCC[current_x][y] = 0;

        }
    } else { //if vertical
        //-> remove the particle
        for (int i = 0; i < particle_length; i++)
        {
            int current_y = y + i;
            if (current_y > M - 1) current_y = y + i - M; //jump to opposite lattice point if exceeding lattice limit
            OCC[x][current_y] = 0;

        }
    }
}

//prints the particles correspondingly for a pyhton visualisation file
static void printParticlesForVisualisation(std::vector<std::vector<int>> rods) {
    //Create files for vertical and horizontal particles
    ofstream vertical("Senkrechte.dat"); 
    ofstream horizontal("Waagerechte.dat");
    for (std::vector<int> element : rods) { //loop throug pool of particles and 
        element[2] ? horizontal << element[0] << " " << element[1] << "\n" : vertical << element[0] << " " << element[1] << "\n"; //write it to corresponding file
    }
}


//comprises a complete GCMC-step, z is activity
static void GCMC(
    bitset<M> (& OCC)[M],
    vector<vector<int>> &rods,
    uniform_int_distribution<int> &unif_distr_pos,
    uniform_real_distribution<double> &unif_distr,
    mt19937 &gen, double z, int& N_v, int& N_h) {
    
    uniform_int_distribution<int> dis1(0,1);

    int N = N_v + N_h; //compute total N
    
    //roll probabilities
    bool insertion = dis1(gen); //true -> insertion, false -> deletion
   
    bool horizontal = dis1(gen); //true -> horizontal, false -> vertical
    
    double r_alpha = unif_distr(gen); //roll for next prob
    //prng for deletion or insertion probability
        //insertion
        if (insertion) { //inserting 
            
            //roll positions selected for insertion
            int r_x = unif_distr_pos(gen);
            int r_y = unif_distr_pos(gen);

            /*
            while (OCC[r_x][r_y] == 1) {
                 r_x = unif_distr_pos(gen);
                 r_y = unif_distr_pos(gen);
            }
           */
            double alpha = min(2 * z * (pow(M, 2) / ((double)(N + 1))), 1.0); //acceptance probability for insertion
           
           if (r_alpha <= alpha) { //try acceptance
               if (add_particle(OCC, { r_x, r_y }, horizontal)) {  //try adding particle, if successful add it to rods
                   rods.push_back({ r_x, r_y, horizontal });
                   horizontal ? N_h++ : N_v++; //increase vertical or horizontal counter
               }
           }else {
               return; //if not successful do nothing
           }
            
        }
        
        else { //deletion
            
            double alpha = min(N/z * 1 / (double)(2 * (pow(M, 2))), 1.0); //acceptance probability for  deletion
            
            if (r_alpha <= alpha) { //try acceptance
                if (N == 0) return;
                //roll random index
                uniform_int_distribution<int> dis(0, N - 1);
                int rand = dis(gen);
                vector<int> props = rods[rand];
                remove_particle(OCC, props); //remove particle from OCC
                rods.erase(rods.begin() + rand); //remove particle froom list
                props[2] ? N_h-- : N_v--; //decrease vertical or horizontal counter
                
            }

        }




}

//starts the simulation with the stated input params
static void simulate(string filename, double z, std::uint64_t thermal_time, std::uint64_t steps_between_measurements, std::uint64_t number_of_measurements, unsigned int seed) {
    mt19937 gen(seed); //mt19937 (mersenne twister) is a very good RNG (has a period of 2^(19937)-1)
    uniform_real_distribution<double> unif_distr(0.0, 1.0); //distribution to roll numbers between 0.0 and 1.0
    uniform_int_distribution<int> unif_distr_pos(0, M - 1); //distribution for rolling x,y values of lattice points


    const int dim_length = M; //length of the lattice box. size -> (MxM)

    bitset<M> OCC[M]; //chosen data structure of OCC-lattice
    //bitset enables faster operations than vector or list. 
    //this code can be optimized even more with using actual bit operations


    init_OCC(OCC, M); //initialzise the lattice

    int N_h = 0; // counter for horizontal particles
    int N_v = 0; //counter for vertical partcles
    vector<vector<int>> rods = {}; //bookkeeping of rods {x_pos, y_pos, orientation}
    

    for (int i = 0; i < thermal_time; i++) //thermalisation
    {
        GCMC(OCC,rods ,unif_distr_pos, unif_distr, gen, z, N_v, N_h);
    }
   
    //write data to file
    ofstream data_file(filename); 
    //write important generation params into file as a header
    data_file << "simulation params: \n" << "z: " << z << "\n" << "thermalisation steps: " << thermal_time << "\n";
    data_file << "no. measurements: " << number_of_measurements << "\n" << "interval: " << steps_between_measurements << "\n";
    data_file << "RNG seed: " << seed << "\n";
   

    data_file << "i" << "\t" << "N" << "\t" << "N_h" << "\t" << "N_v" << "\n";
    for (unsigned int i = 0; i < number_of_measurements; i++) { //total simulation steps: (no. measurements x interval) + therm_time 
        for (unsigned int j = 0; j < steps_between_measurements; j++)
        {
            GCMC(OCC, rods, unif_distr_pos, unif_distr, gen, z, N_v, N_h);
        }
        int N = N_h + N_v;
        data_file << i << "\t" << N << "\t" << N_h << "\t" << N_v << "\n";


    }
   
    printParticlesForVisualisation(rods); //writes the data of the last configuration to 2 seperate files (vertical and horizontal particles)
    //print_OCC(OCC); //can be uncommented to print the last OCC lattice config to the console in Matrix form
   

    data_file.close();




}



int main(int argc, char* argv[])
{
    
    if (argc != 7) { //if the given no. of params are not correct provide a error message
        std::cerr << "Usage:\n";
        std::cerr << "Usage: sim.exe <output_filename> <z> <thermalisation_steps> <interval> <total_measurements> <seed> \n";
        return 1;
    }

    //parsing simulation params from a .bat file
    std::string output_filename = argv[1];
    double z = std::atof(argv[2]);
    
    //parsing large numbers to from a .bat file with std::stoull() (as string)
    std::uint64_t thermalisation_steps = std::stoull(argv[3]);
    std::uint64_t interval = std::stoull(argv[4]);
    std::uint64_t total_measurements = std::stoull(argv[5]);
    unsigned int seed = std::atoi(argv[6]);

    clock_t start = clock();
    simulate(output_filename, z, thermalisation_steps, interval, total_measurements, seed); //start the simulation
    
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "\n";
    printf("execution time for z0 = %.2f: %.3f sec", z, elapsed);
    cout << "\n";


}









