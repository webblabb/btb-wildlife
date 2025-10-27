#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <chrono>
#include <random>
#include <sstream>
#include <algorithm>
#include <math.h>

#define _USE_MATH_DEFINES
#include <cmath>

#include "gsl/gsl_randist.h"

//vestigial functions from WH cattle model that may be of use???
//For splitting a string into a vector based on specified delimiter.
std::vector<std::string> split(const std::string &s,
                               char delim)
{
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elems;
    elems.reserve(100);
    while (std::getline(ss, item, delim))
    {
        if (!item.empty())
        {
            elems.push_back(item);
        }
    }
    return elems;
}

//Removes whitespace around string.
std::string trim(const std::string &s)
{
    auto wsfront=std::find_if_not(s.begin(),s.end(),[](int c){return std::isspace(c);});
    auto wsback=std::find_if_not(s.rbegin(),s.rend(),[](int c){return std::isspace(c);}).base();
    return (wsback<=wsfront ? std::string() : std::string(wsfront,wsback));
}

//Removes the initial byte-order-mark character added by some software (Excel for instace).
void skipBOM(std::ifstream &in)
{
    char test[3] = {0};
    in.read(test, 3);
    if ((unsigned char)test[0] == 0xEF &&
        (unsigned char)test[1] == 0xBB &&
        (unsigned char)test[2] == 0xBF)
    {
        return;
    }
    in.seekg(0);
}

//Get data from standard input sream.
int read_cin(std::vector<std::string>& out)
{
    std::vector<std::string> vals;
    vals.reserve(200);

    std::string val;
    while(std::getline(std::cin, val))
    {
        if(!val.empty())
        {
            vals.push_back(val);
        }
    }
    out.swap(vals);
    return 0;
}

//Read the command line arguments upon execution.
int read_args(int argc, char* argv[], std::vector<std::string>& out)
{
    std::vector<std::string> vals;
    vals.reserve(200);

    std::string val;
    for(int i=1; i<argc; ++i)
    {
        val = argv[i];
        vals.push_back(val);
    }
    out.swap(vals);
    return 0;
}

//
//

struct ReplicateResult
{
    int n_events = 1;
    
    std::vector<double> res_time; //1
    
    std::vector<int> res_n; //2a
    std::vector<int> res_nSuper; //2b
    
    std::vector<int> res_nS; //3a
    std::vector<int> res_nSuperS; //3b
    
    std::vector<int> res_nE1; //4a
    std::vector<int> res_nSuperE1; //4b
    
    std::vector<int> res_nI; //5a
    std::vector<int> res_nSuperI; //5b
    
    std::vector<int> res_Total_inf; //6
    
    std::vector<int> res_month; //7
    std::vector<int> res_quarter; //8
    std::vector<int> res_event;
    
    std::vector<double> res_lambda; //9
    

    ReplicateResult(int n_events)
    {
        res_time.assign(1,-1); //1
        
        res_n.assign(1,-1); //2a
        res_nSuper.assign(1,-1); //2b
        
        res_nS.assign(1,-1); //3a
        res_nSuperS.assign(1,-1); //3b
        
        res_nE1.assign(1,-1); //4a
        res_nSuperE1.assign(1,-1); //4b
        
        res_nI.assign(1,-1); //5a
        res_nSuperI.assign(1,-1); //5b
        
        res_Total_inf.assign(1,-1); //6
        
        res_month.assign(1,-1); //7
        res_quarter.assign(1,-1); //8
        res_event.assign(1,-1);
        
        res_lambda.assign(1,-1); //9
    }

};


//Struct for storing the various parameters that goes into the model.
struct Parameters
{
    //Number of animals to start with in infection classes.
    int S_0 = 200; //Susceptible
    int E1_0 = 0; //Exposed
    int I_0 = 0; //Infectious
    
    //super spreader classes implemented to address disproportionate interactions with farms due to a select few individuals.
    //These animals will only differ from normal individuals of the same class by a scaling factor (phi) on p2.

    int SuperS_0 = 0; //SS Susceptible
    int SuperE1_0 = 0; //SS Exposed
    int SuperI_0 = 0; //SS Infectious

    //herd/county -level parameters
    int herd_size = S_0 + E1_0 + I_0 + SuperS_0 + SuperE1_0 + SuperI_0; // N, true herd size
    
    
    
    //seasonal parameters
    /////////////////////////////////////////
    //mortality
    double carrying_capacity = herd_size; // habitat carrying capacity
    double hunting_rate = .1 / 3.0; // eta_h - not actually a rate...
    double mortality_rate = .1; // eta_n
    double competition_power = 1; // theta
    double competition_scalar = 1; // gamma
    
    
    //seasonal birth pulses
    double max_birth_rate = 1 ; // birth pulse amplitude
    double prop_newborn_superSpreader = .05; //ksi
    double phase_shift = .1; // omega - birth pulse phase shift
    double synchrony = .1; // s - birth pulse synchrony
    double annual_birth_rate = 0; //alpha
    
    
    //////////////////////////////////////////
    
    //Transmision parameters
    double transmission_rate = 25.0; // beta
    double county_habitat_area = 1000.0; //A
    double deer_deer_contact_rate = 0.3; // p1
    double deer_farm_contact_rate_q1 = 0.0; // p2_q1
    double deer_farm_contact_rate_q2 = 0.0; // p2_q2
    double deer_farm_contact_rate_q3 = 0.0; // p2_q3
    double deer_farm_contact_rate_q4 = 0.0; // p2_q4
    double SS_farm_contact_factor = 2.0; //phi

    //sigma_1 is the rate of transition from class E1 to I and is drawn from a gamma distribution
    //with these parameters
    double sigma1_mean = 14.0;// sigma1 //actually mean of lomax distn. that gives the desired gamma mean over all individuals
    double sigma1_rate = 0.5; // sigma1_rate
    double sigma1_scale = 1.0 / sigma1_rate; //Just for use with scipys gamma rv fun that takes scale rather than rate.
    double sigma1_shape = sigma1_mean * sigma1_rate; // sigma1_rate

    //General settings
    int n_years = 5;
    int start_quarter = 4;

    std::map<int, int> start_quarter_to_month = {
        {1, 1}, {2, 4}, {3, 7}, {4, 10}
    };
    
    std::map<int, int> month_to_quarter = {
        {1, 0}, {2, 0}, {3, 0},
        {4, 1}, {5, 1}, {6, 1},
        {7, 2}, {8, 2}, {9, 2},
        {10, 3}, {11, 3}, {12, 3}
    };
    
    int start_month = start_quarter_to_month[start_quarter];
    std::string batchname = "bTB_wl_BS_test";
    
    
    int verbose = 1;
    std::vector<std::string> conf_v;

    //Constructor either reads parameters from standard input (if no argument is given),
    //or from file (argument 1(.
    Parameters(int argc, char* argv[])
    {
        conf_v.reserve(200);
        std::stringstream buffer;
        if(argc == 3) //Config file
        {
            std::ifstream f(argv[2]);
            if(f.is_open())
            {
                buffer << f.rdbuf();
            }
            else
            {
                std::cout << "Failed to read config file \"" << argv[2] << "\"" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            for(int i=2; i<argc; ++i)
            {
                buffer << argv[i] << std::endl;
            }
        }
        while(!buffer.eof())
        { // until end of the stream
            std::string line = "";
            std::stringstream line_ss;
            // First get a whole line from the config file
            std::getline(buffer, line);
            // Put that in a stringstream (required for getline) and use getline again
            // to only read up until the comment character (delimiter).
            line_ss << line;
            std::getline(line_ss, line, '#');
            // If there is whitespace between the value of interest and the comment character
            // this will be read as part of the value. Therefore strip it of whitespace first.
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            if(line.size() != 0)
            {
                conf_v.emplace_back(line);
            }
        }

        size_t conf_length = 30;
        if (conf_v.size()!=conf_length)
        {
            std::cout << "Expected configuration file with " << conf_length << " options, loaded file with "
                      << conf_v.size() << " lines."<<std::endl;
            //std::cout << "first: " << std::stod(conf_v[0]) << std::endl;
            //std::cout << "last: " << std::stod(conf_v[22]) << std::endl;
    
            exit(EXIT_FAILURE);
        }

        //Population parameters
        S_0 = std::stoi(conf_v[0]); //Susceptible
        E1_0 = std::stoi(conf_v[1]); //Exposed
        I_0 = std::stoi(conf_v[2]); //Infectious
        SuperS_0 = std::stoi(conf_v[3]); //Susceptible super spreader
        SuperE1_0 = std::stoi(conf_v[4]); //Exposed super spreader
        SuperI_0 = std::stoi(conf_v[5]); //Infectious super spreader
        
        
        //seasonal parameters
    
        //mortality
        carrying_capacity = std::stoi(conf_v[6]); // K, carrying capacity / start population
        hunting_rate = std::stod(conf_v[7]); // eta_hunt
        mortality_rate = std::stod(conf_v[8]); // eta_nat
        competition_power = std::stod(conf_v[9]); // theta
        competition_scalar = std::stod(conf_v[10]); // gamma
        
        //seasonal birth pulse parameters
        max_birth_rate = std::stod(conf_v[11]); // A, birth pulse amplitude
        prop_newborn_superSpreader = std::stod(conf_v[12]); //ksi, proporton of newborn 'Super Spreaders'
        phase_shift = std::stod(conf_v[13]); // omega, birth pulse phase shift
        synchrony = std::stod(conf_v[14]); // s, birth pulse synchrony
        annual_birth_rate = std::stod(conf_v[15]); // alpha

        //Transmision parameters
        transmission_rate = std::stod(conf_v[16]); // beta
        county_habitat_area = std::stod(conf_v[17]); // Area
        deer_deer_contact_rate = std::stod(conf_v[18]); // p1
        deer_farm_contact_rate_q1 = std::stod(conf_v[19]); // p2_q1
        deer_farm_contact_rate_q2 = std::stod(conf_v[20]); // p2_q2
        deer_farm_contact_rate_q3 = std::stod(conf_v[21]); // p2_q3
        deer_farm_contact_rate_q4 = std::stod(conf_v[22]); // p2_q4
        SS_farm_contact_factor = std::stod(conf_v[23]); // phi
        
        //sigma_1 is the rate of transition from class E1 to I and is drawn from a gamma distribution
        //with these parameters
        sigma1_mean = std::stod(conf_v[24]);// sigma1
        sigma1_rate = std::stod(conf_v[25]); // sigma1_rate
        sigma1_scale = 1.0 / sigma1_rate; //Just for use with scipys gamma rv fun that takes scale rather than rate.
        sigma1_shape = sigma1_mean * sigma1_rate; // sigma1_rate


        //General settings
        n_years = std::stoi(conf_v[26]);
        start_quarter = std::stoi(conf_v[27]);
        start_month = start_quarter_to_month.at(start_quarter);
        verbose = std::stod(conf_v[28]);
        batchname = conf_v[29];
    }
};

//For initiating and seeding gsl random number generator.
class Gsl_rng_wrapper
{
    gsl_rng* r;
    public:
        Gsl_rng_wrapper()
        {
            std::random_device rng_dev; //To generate a safe seed.
            long seed = time(NULL)*rng_dev();
            const gsl_rng_type* rng_type = gsl_rng_default;
            r = gsl_rng_alloc(rng_type);
            gsl_rng_set(r, seed);
        }
        ~Gsl_rng_wrapper() { gsl_rng_free(r); }
        gsl_rng* get_r() { return r; }
};

//Uniform RV using gsl.
double draw_uniform_gsl(double lo, double hi)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_flat(r, lo, hi);
}

//Binomial RV using gsl.
int draw_binom_gsl(int N, double prob)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_binomial(r, prob, N);
}

//Gamma RV using gsl.
double draw_gamma_gsl(double shape, double scale)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_gamma(r, shape, scale);
}

//Poisson RV using gsl
int draw_poisson_gsl(double rate)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_poisson(r, rate);
}

/*/////////////////////////////////
//////////////////////////////////
////Full Wildlife method v.1.0///
////////////////////////////////
/////////////////////////////*/


//continuous time model CTMC
ReplicateResult btb_wl_model_CTMC(const Parameters& p)
{
    int verbose = p.verbose;
    if(verbose > 3){std::cout << "begin function body" << std::endl;}
    
    
    if(verbose > 1){std::cout << "begin parameter assignment" << std::endl;}
    
        //create vectors for seasonal parameters
        std::vector<double> farm_contact_rate_v = {
            p.deer_farm_contact_rate_q1,
            p.deer_farm_contact_rate_q2,
            p.deer_farm_contact_rate_q3,
            p.deer_farm_contact_rate_q4
        };


        int n_tsteps = p.n_years * 12; //total number of timesteps
        double tstep_frac = double(p.n_years) / double(n_tsteps);
    
    
    if(verbose > 1){std::cout << "begin res object creation" << std::endl;}
    
        ReplicateResult res(1);
    if(verbose > 1){std::cout << "end res object creation " << res.res_event[0] << res.res_lambda[0] << res.res_lambda.size() << std::endl;}
    
    if(verbose > 1){std::cout << "initializing state variables" << std::endl;}
    
        int nS = p.S_0; //3a
        int nSuperS = p.SuperS_0; //3b
    
        int nE1 = p.E1_0; //4a
        int nSuperE1 = p.SuperE1_0; //4b
    
        int nI = p.I_0; //5a
        int nSuperI = p.SuperI_0; //5b
    
        int N = nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI; //2a
        int NSuper = nSuperS + nSuperE1 + nSuperI;  //2b
    
        double N_d = nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI; //N defined as a double for comupting host frequency
        int event = 0;
        int tstep_0 = 0;
        double tot_rate = 0;
    
    
    if(verbose > 1){std::cout << "assigning state variable values" << std::endl;}
    
        res.res_time[tstep_0] = tstep_0; //1
    
        res.res_n[tstep_0] = N; //2a
        res.res_nSuper[tstep_0] = NSuper; //2b
    
        res.res_nS[tstep_0] = nS; //3a
        res.res_nSuperS[tstep_0] = nSuperS; //3b
    
        res.res_nE1[tstep_0] = nE1; //4a
        res.res_nSuperE1[tstep_0] = nSuperE1; //4b
    
        res.res_nI[tstep_0] = nI; //5a
        res.res_nSuperI[tstep_0] = nSuperI; //5b
    
        res.res_event[tstep_0] = 0; //6
    
        res.res_Total_inf[tstep_0] = nE1 + nSuperE1 + nI + nSuperI;
        res.res_month[tstep_0] = p.start_month; //7
        res.res_quarter[tstep_0] = p.month_to_quarter.at(p.start_month) + 1; //8
   
        res.res_lambda[tstep_0] = 0; //9
    
    if(verbose > 2){std::cout << "task complete" << std::endl;}
    
    if(verbose > 0){std::cout << "beginning loop over all timesteps" << std::endl;}
    
    //start loop over all timesteps
    double time = 0;
    while(time < n_tsteps)
    {
        if(verbose > 3){std::cout << "start timestep " << time << std::endl;}
        
        
        //assign month, quarter, and then seasonal variables
        if(verbose > 1){std::cout << "task: assign seasonal variables" << std::endl;}
            int current_month = int( floor(time) + p.start_month - 1) % 12 + 1;
        
        if(verbose > 5){std::cout << "month: " << current_month << std::endl;}
        
            int current_quarter_idx = p.month_to_quarter.at(current_month);
        if(verbose > 5){std::cout << "month index: " << current_quarter_idx << std::endl;}
        
            int current_quarter = p.month_to_quarter.at(current_month) + 1;
        if(verbose > 5){std::cout << "quarter: " << current_quarter << std::endl;}
            N_d = double(N);
        
        int S_change = 0;
        int S_SS_change = 0;
        int E1_change = 0;
        int E1_SS_change = 0;
        int I_change = 0;
        int I_SS_change = 0;
        //should only take values of -1, 0, 1 since only 1 event can occur per unit time
        
            //compute rates for all events
        
            //hunting
            double hunting_rate = 0.0;
            if(current_quarter == 4)
            {
                hunting_rate = p.hunting_rate;
            }
        
            //natural mortality
            double nat_mortality_rate = (N_d)*(p.mortality_rate + p.competition_scalar * pow((N_d / double(p.carrying_capacity)), p.competition_power));
            
        
            //birth
            double BR = (N_d)*p.max_birth_rate*exp( -p.synchrony*pow( cos( (M_PI/12) * (time + double(p.start_month-1) - p.phase_shift) ) , 2.0 ));
    
            //transmission, separated by mechanism
            
            //deer-deer transmission (density dependence)
            double interspecific_trans_rate = p.deer_deer_contact_rate * p.transmission_rate * (double(nS+nSuperS)) * (double(nI+nSuperI) / (p.county_habitat_area));
        
            //transmission from infected farm contact
            double farm_trans_rate = farm_contact_rate_v[current_quarter_idx] * p.transmission_rate * (nS+p.SS_farm_contact_factor*nSuperS);
        
            //transition
            double transition_rate = (nE1 + nSuperE1) / (draw_gamma_gsl(p.sigma1_shape, p.sigma1_scale));
            if(p.sigma1_shape == 0 || p.sigma1_scale == 0){transition_rate = 0;}
            //determine total event rate
            tot_rate = BR + hunting_rate + nat_mortality_rate + interspecific_trans_rate + farm_trans_rate + transition_rate;
            //ordered birth, natural mortality, hunt mortality, transmission, transition
            
        if(verbose > 1){
            std::cout << "printing rates: " << std::endl;
            std::cout << "hunt rate: " << hunting_rate << std::endl;
            std::cout << "mort. rate: " << nat_mortality_rate << std::endl;
            std::cout << "Birth rate: " << BR << std::endl;
            std::cout << "transmission rate: " << interspecific_trans_rate << std::endl;
            std::cout << "spill. transmission rate: " << farm_trans_rate << std::endl;
            std::cout << "transition rate: " << transition_rate << std::endl;
            
        }
        
        //determine time to next event
        double scale_frac = 1 / tot_rate;
        
        if( (N == 0 && NSuper == 0)  ){
            scale_frac = 1;
            tot_rate = 0;
            BR = 0;
            hunting_rate = 0;
            nat_mortality_rate = 0;
            interspecific_trans_rate = 0;
            farm_trans_rate = 0;
            transition_rate = 0;
        }
        
        double tau = -(std::log(draw_uniform_gsl(0.0, 1.0))) * scale_frac;
        if(verbose > 5){std::cout << "time: " << time << " + " << tau << std::endl;}
        
        //control loop for testing single events
        bool control = false;
        if(tau > 1.0 and verbose == -1 ){
            tau = 1.0/4.0;
            control = true;
        }
            time += tau; //increment time by the stochastic time to next event
            current_month = int( floor(time) + p.start_month - 1) % 12 + 1;
            
            //determine event type
            double U1; double U2;
            double event_cut = 0.0;
            U1 = draw_uniform_gsl(0.0, 1.0);
        
            /*/////////
               birth
            /////////*/
        if( double(nS + nSuperS + nE1 + nSuperE1 + nI + nSuperI) / N_d < 1){
            std::cout << "Invalid population proportion" << double(nS + nSuperS + nE1 + nSuperE1 + nI + nSuperI) / N_d << std::endl;
            
        }
        
        if(control)
        {
            event = 0;
        }
        
        else if(U1 < (BR)*scale_frac)
        {
            event_cut = (BR)*scale_frac;
            if(verbose > 6){std::cout << "computing birth..." << std::endl;}
            U2 = draw_uniform_gsl(0.0, 1.0);
            
            if( U2 < p.prop_newborn_superSpreader ){ event = 2; S_SS_change = 1;}//S_SSbirth
            else if( U2 < 1 ){event = 1; S_change = 1; }//S birth
            else{event = 0;}
        }
        
            /*///////
              death
            ///////*/
        
            //nat mortality
        else if(U1 < (BR + nat_mortality_rate)*scale_frac)
        {
            event_cut = (BR + nat_mortality_rate)*scale_frac;
            //determine class for mortality
            if(verbose > 6){std::cout << "computing nat mortality..." << std::endl;}
            U2 = draw_uniform_gsl(0.0, 1.0);
            
            if( U2 < (double(nS) / N_d) ){event = 3; S_change = -1; }//S
            else if( U2 < (double(nS + nE1) / N_d) ){event = 4; E1_change = -1; }//E1
            else if( U2 < (double(nS + nE1 + nI) / N_d) ){event = 5; I_change = -1; }//I
            else if( U2 < (double(nS + nE1 + nI + nSuperS) / N_d) ){event = 6; S_SS_change = -1; }//S_SS
            else if( U2 < (double(nS + nE1 + nI + nSuperS + nSuperE1) / N_d) ){event = 7; E1_SS_change = -1; }//E1_SS
            else if( U2 < (double(nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI) / N_d) ){event = 8; I_SS_change = -1; }//I_SS
            else{event = 0;}
        }
        
            //hunt mortality
        else if(U1 < (BR + nat_mortality_rate + hunting_rate)*scale_frac)
        {
            event_cut = (BR + nat_mortality_rate + hunting_rate)*scale_frac;
            //determine hunted class
            if(verbose > 6){std::cout << "computing hunt mortality..." << std::endl;}
            U2 = draw_uniform_gsl(0.0, 1.0);
            
            if( U2 < (double(nS) / N_d) ){event = 9; S_change = -1; }//S
            else if( U2 < (double(nS + nE1) / N_d) ){event = 10; E1_change = -1; }//E1
            else if( U2 < (double(nS + nE1 + nI) / N_d) ){event = 11; I_change = -1; }//I
            else if( U2 < (double(nS + nE1 + nI + nSuperS) / N_d) ){event = 12; S_SS_change = -1; }//S_SS
            else if( U2 < (double(nS + nE1 + nI + nSuperS + nSuperE1) / N_d) ){event = 13; E1_SS_change = -1; }//E1_SS
            else if( U2 < (double(nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI) / N_d) ){event = 14; I_SS_change = -1;}//I_SS
            else{event = 0;}
        }
        
        
            /*/////////////////
                Transmission
            /////////////////*/
        
        //deer-deer contact transmission
        
        else if(U1 < (BR  + nat_mortality_rate + hunting_rate + interspecific_trans_rate)*scale_frac)
        {
            event_cut = (BR  + nat_mortality_rate + hunting_rate + interspecific_trans_rate)*scale_frac;
            //divide between normal and SS susceptible individuals
            if(verbose > 6){std::cout << "computing i transmission..." << std::endl;}
            U2 = draw_uniform_gsl(0.0, 1.0);
            
            if( U2 < (double(nS)/ double(nS+nSuperS)) ){event = 15; S_change = -1; E1_change = 1;} //S becomes infected
            else if( U2 < ( double(nS+nSuperS)/double(nS+nSuperS) ) ){event = 16; S_SS_change = -1; E1_SS_change = 1;} //S_SS individual becomes infected
            else{event = 0;}
        }
        
        
        //deer-farm/fomite contact transmission
        
        else if(U1 < (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate + farm_trans_rate)*scale_frac)
        {
            event_cut = (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate + farm_trans_rate)*scale_frac;
            //divide between normal and SS individuals
            if(verbose > 6){std::cout << "computing s transmission..." << std::endl;}
            U2 = draw_uniform_gsl(0.0, 1.0);
            if( U2 < (double(nS) / double(nS+p.SS_farm_contact_factor*nSuperS)) ){event = 17; S_change = -1; E1_change = 1;} //S becomes infected
            else if( U2 < double(nS+p.SS_farm_contact_factor*nSuperS) / double(nS+p.SS_farm_contact_factor*nSuperS) ){event = 18; S_SS_change = -1; E1_SS_change = 1;} //S_SS individual becomes infected
            else{event = 0;}
        }
       
                        
        /*//////////////
            Transition
        //////////////*/
        
        else if(U1 < (BR + hunting_rate + nat_mortality_rate + interspecific_trans_rate + farm_trans_rate + transition_rate)*scale_frac){
            event_cut = (BR + hunting_rate + nat_mortality_rate + interspecific_trans_rate + farm_trans_rate + transition_rate)*scale_frac;
            //divide between normal or SS exposed individuals
            if(verbose > 6){std::cout << "computing transition..." << std::endl;}
            U2 = draw_uniform_gsl(0.0, 1.0);
            
            if( U2 < (double(nE1) / double(nE1+nSuperE1)) ){event = 19; E1_change = -1; I_change = 1;} //E1 becomes infectious
            else if( U2 < (double(nE1+nSuperE1) / double(nE1+nSuperE1)) ){event = 20; E1_SS_change = -1; I_SS_change = 1;} //E1_SS individual becomes infectious
            else{event = 0;}
        }

        else{
            if(verbose > 6){std::cout << "No event..." << std::endl;}
            event = 0;
        }
        
        
        if(verbose > 6){
            std::cout << "Upper event cutoffs: " << std::endl;
            std::cout << "U1: " << U1 << std::endl;
            std::cout << "Birth rate: " << (BR)*scale_frac << std::endl;
            std::cout << "mort. rate: " << (BR + nat_mortality_rate)*scale_frac << std::endl;
            std::cout << "hunt rate: " << (BR + nat_mortality_rate + hunting_rate)*scale_frac << std::endl;
            std::cout << "transmission rate: " << (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate)*scale_frac << std::endl;
            std::cout << "spill. transmission rate: " << (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate + farm_trans_rate)*scale_frac << std::endl;
            std::cout << "transition rate: " << (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate + farm_trans_rate + transition_rate)*scale_frac << std::endl;
        }
        
       //Update classes
        
        nS += S_change;
        nSuperS += S_SS_change;
        nE1 += E1_change;
        nSuperE1 += E1_SS_change;
        nI += I_change;
        nSuperI += I_SS_change;
        N = nS + nSuperS + nE1 + nSuperE1 + nI + nSuperI;
        NSuper = nSuperS + nSuperE1 + nSuperI;
        
        if(nS < 0 || nSuperS < 0 || nE1 < 0 || nSuperE1 < 0 || nI < 0 || nSuperI < 0 || N < 0 || NSuper < 0){
            std::cout << "Negative animals: " << std::endl;
            std::cout << "Event: " << event << std::endl;
            std::cout << "nS: " << nS << " change: " << S_change << std::endl;
            std::cout << "nSuperS: " << nSuperS << " change: " << S_SS_change << std::endl;
            std::cout << "nE1: " << nE1 << " change: " << E1_change << std::endl;
            std::cout << "nSuperE1: " << nSuperE1 << " change: " << E1_SS_change << std::endl;
            std::cout << "nI: " << nI << " change: " << I_change << std::endl;
            std::cout << "nSuperI: " << nSuperI << " change: " << I_SS_change << std::endl;
            std::cout << "N: " << N << std::endl;
            std::cout << "NSuper: " << NSuper << std::endl;
            std::cout << "U1: " << U1 << std::endl;
            std::cout << "U2: " << U2 << std::endl;
            std::cout << "N " << N_d << std::endl;
            std::cout << double(nS + nE1 + nI + nSuperS + nSuperE1) / N_d << std::endl;
            break;
        }
        
        if(time > n_tsteps){time = n_tsteps;}//
        
        //apply changes to res structure
        res.res_time.push_back(time); //1
    
        res.res_n.push_back(N); //2a
        res.res_nSuper.push_back(NSuper); //2b
    
        res.res_nS.push_back(nS); //3a
        res.res_nSuperS.push_back(nSuperS); //3b
    
        res.res_nE1.push_back(nE1); //4a
        res.res_nSuperE1.push_back(nSuperE1); //4b
    
        res.res_nI.push_back(nI); //5a
        res.res_nSuperI.push_back(nSuperI); //5b
    
        res.res_month.push_back(current_month); //6
        res.res_quarter.push_back( p.month_to_quarter.at(current_month) + 1); //7
        res.res_event.push_back(event);
        
        res.res_Total_inf.push_back(nE1 + nSuperE1 + nI + nSuperI); //8
    
        res.res_lambda.push_back(tot_rate); //9
        
        //std::cout << "end res object creation " << res.res_event[res.res_lambda.size()-1] << " " << res.res_time[res.res_lambda.size()-1] << " " << res.res_lambda.size() << std::endl;
        int brk = verbose;
        if( brk == 606 ){
            if(event > 0 ){
                std::cout << "0 animals detected: " << std::endl;
                std::cout << "Event: " << event << std::endl;
                std::cout << "nS: " << nS << " change: " << S_change << std::endl;
                std::cout << "nSuperS: " << nSuperS << " change: " << S_SS_change << std::endl;
                std::cout << "nE1: " << nE1 << " change: " << E1_change << std::endl;
                std::cout << "nSuperE1: " << nSuperE1 << " change: " << E1_SS_change << std::endl;
                std::cout << "nI: " << nI << " change: " << I_change << std::endl;
                std::cout << "nSuperI: " << nSuperI << " change: " << I_SS_change << std::endl;
                std::cout << "NSuper: " << NSuper << std::endl;
                std::cout << "N: " << N << std::endl;
                std::cout << "Scale fraction: " << scale_frac << std::endl;
                std::cout << "Inv. scale fraction: " << 1.0/scale_frac << std::endl;
                std::cout << "Event sum: " << BR + hunting_rate + nat_mortality_rate + interspecific_trans_rate + farm_trans_rate + transition_rate << std::endl;
                std::cout << "event cutoff: " << event_cut << std::endl;
                std::cout << "U1: " << U1 << std::endl;
                std::cout << "U2: " << U2 << std::endl;
                std::cout << double(nS + nE1 + nI + nSuperS + nSuperE1) / N_d << std::endl;
                
                std::cout << "Upper event cutoffs: " << std::endl;
                std::cout << "Birth rate: " << (BR)*scale_frac << std::endl;
                std::cout << "mort. rate: " << (BR + nat_mortality_rate)*scale_frac << std::endl;
                std::cout << "hunt rate: " << (BR + nat_mortality_rate + hunting_rate)*scale_frac << std::endl;
                std::cout << "transmission rate: " << (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate)*scale_frac << std::endl;
                std::cout << "spill. transmission rate: " << (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate + farm_trans_rate)*scale_frac << std::endl;
                std::cout << "transition rate: " << (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate + farm_trans_rate + transition_rate)*scale_frac << std::endl;
            }
        }
        if( N == 0 && NSuper == 0 ){brk = 606;}
        if( N == 0 && NSuper == 0 && brk > 8008){
            time = n_tsteps;
            res.res_time.push_back(time); //1
        
            res.res_n.push_back(N); //2a
            res.res_nSuper.push_back(NSuper); //2b
        
            res.res_nS.push_back(nS); //3a
            res.res_nSuperS.push_back(nSuperS); //3b
        
            res.res_nE1.push_back(nE1); //4a
            res.res_nSuperE1.push_back(nSuperE1); //4b
        
            res.res_nI.push_back(nI); //5a
            res.res_nSuperI.push_back(nSuperI); //5b
        
            res.res_month.push_back(current_month); //6
            res.res_quarter.push_back( p.month_to_quarter.at(current_month) + 1); //7
            res.res_event.push_back(event);
            
            res.res_Total_inf.push_back(nE1 + nSuperE1 + nI + nSuperI); //8
        
            res.res_lambda.push_back(tot_rate); //9
            break;
        }
        
        
    }
    
    return res;
    
}

         
int main(int argc, char* argv[]){
    /*
     If you run this from the command line, you should provide two arguments - the
     number of replicates and a config file (argc==3). In that case the results will
     be output to a file. The alternative is to provide all the parameters in a single
     long string instead of the config file. This is done by the R-wrapper. In that
     case the results will be written to the standard output stream which is caught
     from within R and saved to a table.
     */
    
    
    if(argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <n replicates> <config file / string of parameter values>" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    int n_reps = std::stoi(argv[1]);
    Parameters p(argc, argv);
    
    std::ostream* out_stream;
    std::string out_fn = p.batchname + "_results.txt";
    std::ofstream out_fs;
    if(argc == 3) //Only use an output file if running with config file. Otherwise output res to stdout.
    {
        out_fs.open(out_fn);
        if(!out_fs.is_open())
        {
            std::cout << "Failed to open output file " << out_fn << "." << std::endl;
            exit(EXIT_FAILURE);
        }
        out_stream = &out_fs;
    }
    else
    {
        out_stream = &std::cout;
    }
    
    /*Write header to output file. Note that The first line for each replicate are the starting parameters when t=0. The following lines are the results at the end of the corresponding time steps starting at 1. Therefore there will be 61 lines per replicate if running 60 time steps etc.*/
    std::stringstream ss;
    
    //std::cout << "checkpoint 1" << std::endl;
    
    
    //Troubleshooting parameter values from wrapper code
    if(p.verbose > 0){
        std::cout << "printing parameters: " << std::endl;
        std::cout << "susc. " << p.S_0 << std::endl;//0
        std::cout << "E1 " << p.E1_0 << std::endl;//1
        std::cout << "I " << p.I_0 << std::endl;//2
        std::cout << "Super S " << p.SuperS_0 << std::endl;//3
        std::cout << "Super E1 " << p.SuperE1_0 << std::endl;//4
        std::cout << "Super I " << p.SuperI_0 << std::endl;//5
        
        std::cout << "K " << p.carrying_capacity << std::endl;//6
        std::cout << "eta_h " << p.hunting_rate << std::endl;//6
        std::cout << "eta_nat " << p.mortality_rate << std::endl;//6
        std::cout << "theta " << p.competition_power << std::endl;//6
        std::cout << "gamma " << p.competition_scalar << std::endl;//6
        
        std::cout << "BR_max " << p.max_birth_rate << std::endl;//7
        std::cout << "ksi " << p.prop_newborn_superSpreader << std::endl;//16
        std::cout << "omega " << p.phase_shift << std::endl;//17
        std::cout << "sync. " << p.synchrony << std::endl;//7
        std::cout << "alpha " << p.annual_birth_rate << std::endl;//16
    
        std::cout << "beta " << p.transmission_rate << std::endl;//22
        std::cout << "area " << p.county_habitat_area << std::endl;//22
        std::cout << "p1 " << p.deer_deer_contact_rate << std::endl;//22
        std::cout << "DF1 " << p.deer_farm_contact_rate_q1 << std::endl;//19
        std::cout << "DF2 " << p.deer_farm_contact_rate_q2 << std::endl;//20
        std::cout << "DF3 " << p.deer_farm_contact_rate_q3 << std::endl;//21
        std::cout << "DF4 " << p.deer_farm_contact_rate_q4 << std::endl;//22
        std::cout << "SS fact. " << p.SS_farm_contact_factor << std::endl;//22
        std::cout << "mean " << p.sigma1_mean << std::endl;//23
        std::cout << "rate " << p.sigma1_rate << std::endl;//24
        std::cout << "yrs. " << p.n_years << std::endl;//25
        std::cout << "Start q. " << p.start_quarter << std::endl;//26
        std::cout << "name " << p.batchname << std::endl;//27
        std::cout << "reps " << n_reps << std::endl;
    }
    
    

    //print out only lambda and month values
    if(p.verbose == -3){
        
        //simulation peoperties
        ss
        << "Month;"
        << "Lambda"
        << std::endl;
        *out_stream << ss.rdbuf();
        
        for(int rep_i=0; rep_i<n_reps; ++rep_i){

            ReplicateResult res = btb_wl_model_CTMC(p);
            int n_tsteps = res.res_event.size();
            std::stringstream rep_ss;
        
        
            for(int i=0; i<n_tsteps; ++i){
            
                rep_ss << res.res_month[i] << ";"
                << res.res_lambda[i]
                << std::endl;
                
            }
            *out_stream << rep_ss.rdbuf();
        }
    }
    else{
        
        //simulation peoperties
        ss  << "rep;" << "tstep;" << "time;"
    
        //Population properties
        << "N;" << "N Super;"
        << "S;" << "SS_S;"
        << "E1;" << "SS_E1;"
        << "I;" << "SS_I;"
    
        << "Total Inf;"
        
        //misc. properties
        << "Event;"
        << "Quarter;"
        << "Month;"
        << "Lambda"
        << std::endl;
        *out_stream << ss.rdbuf();
        
        for(int rep_i=0; rep_i<n_reps; ++rep_i){

            ReplicateResult res = btb_wl_model_CTMC(p);
            int n_tsteps = res.res_event.size();
            std::stringstream rep_ss;
        
            for(int i=0; i<n_tsteps; ++i){
                
                rep_ss << rep_i+1 << ";" //0*
                << i << ";" //1
                    
                << res.res_time[i] << ";"
                
                << res.res_n[i] << ";" //2a
                << res.res_nSuper[i] << ";" //2b
                    
                << res.res_nS[i] << ";" //3a
                << res.res_nSuperS[i] << ";" //3b
                    
                << res.res_nE1[i] << ";"//4a
                << res.res_nSuperE1[i] << ";"//4b
                    
                << res.res_nI[i] << ";" //5a
                << res.res_nSuperI[i] << ";" //5b
                    
                << res.res_Total_inf[i] << ";"
                    
                << res.res_event[i] << ";"
                << res.res_quarter[i] << ";"
                << res.res_month[i] << ";"
                << res.res_lambda[i]
                << std::endl;
            }
            *out_stream << rep_ss.rdbuf();
        }
    }
    return 0;
}


