#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <chrono>
#include <random>
#include <sstream>
#include <algorithm>

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
    int n_tsteps; //simulation length
    std::vector<int> res_tstep; //1
    
    std::vector<int> res_n; //2a
    std::vector<int> res_nSuper; //2b
    
    std::vector<int> res_nS; //3a
    std::vector<int> res_nSuperS; //3b
    
    std::vector<int> res_nE1; //4a
    std::vector<int> res_nSuperE1; //4b
    
    std::vector<int> res_nI; //5a
    std::vector<int> res_nSuperI; //5b
    
    std::vector<int> res_S_hunt; //***a
    std::vector<int> res_SuperS_hunt; //***b
    
    std::vector<int> res_E1_hunt; //6a
    std::vector<int> res_SuperE1_hunt; //6b
    
    std::vector<int> res_I_hunt; //7a
    std::vector<int> res_SuperI_hunt; //7b
    
    //susceptible birth/death
    std::vector<int> res_nS_birth; //8a
    std::vector<int> res_nSuperS_birth; //8b
    
    std::vector<int> res_nS_death; //9a
    std::vector<int> res_nSuperS_death; //9b
    
    //E1 birth/death
    std::vector<int> res_nE1_death; //10a
    std::vector<int> res_nSuperE1_death; //10b
    
    //I birth/death
    std::vector<int> res_nI_death; //11a
    std::vector<int> res_nSuperI_death; //11b
    
    // number of transitions
    std::vector<int> res_trans_to_E1; //12a
    std::vector<int> res_trans_to_SuperE1; //12b
    std::vector<int> res_trans_to_I; //13a
    std::vector<int> res_trans_to_SuperI; //13b
    
    std::vector<int> res_month; //14
    std::vector<int> res_quarter; //15
    
    std::vector<int> res_Total_inf;
    std::vector<int> res_Total_hunt;
    
    std::vector<double> res_inf_hunt; //16
    
    std::vector<double> res_FO;
    std::vector<double> res_FO_time;
    std::vector<double> res_Hunt_prev;
    
    

    ReplicateResult(int n_tsteps) : n_tsteps(n_tsteps)
    {
        res_tstep.resize(n_tsteps, -1); //1
        
        res_n.resize(n_tsteps, -1); //2a
        res_nSuper.resize(n_tsteps, -1); //2b
        
        res_nS.resize(n_tsteps, -1); //3a
        res_nSuperS.resize(n_tsteps, -1); //3b
        
        res_nE1.resize(n_tsteps, -1); //4a
        res_nSuperE1.resize(n_tsteps, -1); //4b
        
        res_nI.resize(n_tsteps, -1); //5a
        res_nSuperI.resize(n_tsteps, -1); //5b
        
        res_S_hunt.resize(n_tsteps, -1); //***a
        res_SuperS_hunt.resize(n_tsteps, -1); //***b
        
        res_E1_hunt.resize(n_tsteps, -1); //6a
        res_SuperE1_hunt.resize(n_tsteps, -1); //6b
        
        res_I_hunt.resize(n_tsteps, -1); //7a
        res_SuperI_hunt.resize(n_tsteps, -1); //7b
        
        res_month.resize(n_tsteps, -1); //8
        res_quarter.resize(n_tsteps, -1); //9
        
        res_Total_inf.resize(n_tsteps, -1);
        res_Total_hunt.resize(n_tsteps, -1);
        
        res_inf_hunt.resize(n_tsteps, -1.0); //10
        
        res_nS_birth.resize(n_tsteps, -1); //11a
        res_nSuperS_birth.resize(n_tsteps, -1); //11b
        
        res_nS_death.resize(n_tsteps, -1); //12a
        res_nSuperS_death.resize(n_tsteps, -1); //12b
        
        res_nE1_death.resize(n_tsteps, -1); //13a
        res_nSuperE1_death.resize(n_tsteps, -1); //13b
        
        res_nI_death.resize(n_tsteps, -1); //14a
        res_nSuperI_death.resize(n_tsteps, -1); //14b
        
        res_trans_to_E1.resize(n_tsteps, -1); //15a
        res_trans_to_SuperE1.resize(n_tsteps, -1); //15b
        
        res_trans_to_I.resize(n_tsteps, -1); //16a
        res_trans_to_SuperI.resize(n_tsteps, -1); //16b
        
        res_quarter.resize(n_tsteps, -1);
        
        res_FO.resize(n_tsteps, -1);
        res_FO_time.resize(n_tsteps, -1);
        res_Hunt_prev.resize(n_tsteps, -1);
    
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

    //monthly lambda values
    double L1 = 0.0;
    double L2 = 0.0;
    double L3 = 0.0;
    double L4 = 0.0;
    double L5 = 0.0;
    double L6 = 0.0;
    double L7 = 0.0;
    double L8 = 0.0;
    double L9 = 0.0;
    double L10 = 0.0;
    double L11 = 0.0;
    double L12 = 0.0;
    
    //General settings
    int n_years = 5;
    int start_quarter = 4;
    int integrate_type = 0;
    

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

        size_t conf_length = 43;
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
        
        //monthly lambda values
        L1 = std::stod(conf_v[26]);
        L2 = std::stod(conf_v[27]);
        L3 = std::stod(conf_v[28]);
        L4 = std::stod(conf_v[29]);
        L5 = std::stod(conf_v[30]);
        L6 = std::stod(conf_v[31]);
        L7 = std::stod(conf_v[32]);
        L8 = std::stod(conf_v[33]);
        L9 = std::stod(conf_v[34]);
        L10 = std::stod(conf_v[35]);
        L11 = std::stod(conf_v[36]);
        L12 = std::stod(conf_v[37]);

        //General settings
        integrate_type = std::stoi(conf_v[38]);
        n_years = std::stoi(conf_v[39]);
        start_quarter = std::stoi(conf_v[40]);
        start_month = start_quarter_to_month.at(start_quarter);
        verbose = std::stod(conf_v[41]);
        batchname = conf_v[42];
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

//discrete time model DTMC
ReplicateResult btb_wl_model_DTMC(const Parameters& p)
{
    //preset verbose level until added to parameter input
    //int verbose = 5;
    int verbose = p.verbose;
    
    
    if(p.verbose > 3){std::cout << "begin function body" << std::endl;}
    
    
    if(p.verbose > 3){std::cout << "begin parameter assignment" << std::endl;}
    
    //create vectors for seasonal parameters
    std::vector<double> farm_contact_rate_v = {
        p.deer_farm_contact_rate_q1,
        p.deer_farm_contact_rate_q2,
        p.deer_farm_contact_rate_q3,
        p.deer_farm_contact_rate_q4
    };
    
    //create vectors for seasonal parameters
    std::vector<double> lambda_month = {
        p.L1, p.L2, p.L3, p.L4, p.L5, p.L6,
        p.L7, p.L8, p.L9, p.L10, p.L11, p.L12
    };

    int n_tsteps = p.n_years * 12; //total number of timesteps
    double tstep_frac = double(p.n_years) / double(n_tsteps);
    
    if(p.verbose > 3){std::cout << "task complete" << std::endl;}
    
    
    if(p.verbose > 2){std::cout << "begin res object creation" << std::endl;}
    
        ReplicateResult res(n_tsteps + 1);
    
    if(p.verbose > 2){std::cout << "task complete" << std::endl;}
    
    
    if(p.verbose > 2){std::cout << "initializing state variables" << std::endl;}
    
    int nS = p.S_0; //3a
    int nSuperS = p.SuperS_0; //3b
    
    int nE1 = p.E1_0; //4a
    int nSuperE1 = p.SuperE1_0; //4b
    
    int nI = p.I_0; //5a
    int nSuperI = p.SuperI_0; //5b
    
    int N = nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI; //2a
    int NSuper = nSuperS + nSuperE1 + nSuperI;  //2b
    
    int nS_deaths_hunt = 0; //***a
    int nSuperS_deaths_hunt = 0; //***b
    
    int nE1_deaths_hunt = 0; //6a
    int nSuperE1_deaths_hunt = 0; //6b
    
    int nI_deaths_hunt = 0; //7a
    int nSuperI_deaths_hunt = 0; //7b
    
    double N_d = nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI; //N defined as a double for comupting host frequency
    
    int nS_births = 0; //8
    int nSuperS_births = 0; //8b
    
    
    int nS_deaths = 0; //9a
    int nSuperS_deaths = 0; //9b
    
    int nE1_deaths = 0; //10a
    int nSuperE1_deaths = 0; //10b
    
    int nI_deaths = 0; //11a
    int nSuperI_deaths = 0; //11b
    
    int n_trans_to_E1 = 0; //12a
    int n_trans_to_SuperE1 = 0; //12b
    
    int n_trans_to_I = 0; //13a
    int n_trans_to_SuperI = 0; //13b
    
    //14, 15 i.e. month and quarter are calculated within the loop
    
    double BR = 0; //16
    int tstep_0 = 0;
    int FO = 0;
    int FO_time = 0;
    
    if(p.verbose > 2){std::cout << "task complete" << std::endl;}
    
    
    if(p.verbose > 2){std::cout << "assigning state variable values" << std::endl;}
    
        res.res_tstep[tstep_0] = tstep_0; //1
    
        res.res_n[tstep_0] = N; //2a
        res.res_nSuper[tstep_0] = NSuper; //2b
    
        res.res_nS[tstep_0] = nS; //3a
        res.res_nSuperS[tstep_0] = nSuperS; //3b
    
        res.res_nE1[tstep_0] = nE1; //4a
        res.res_nSuperE1[tstep_0] = nSuperE1; //4b
    
        res.res_nI[tstep_0] = nI; //5a
        res.res_nSuperI[tstep_0] = nSuperI; //5b
    
        res.res_S_hunt[tstep_0] = nS_deaths_hunt; //***a
        res.res_SuperS_hunt[tstep_0] = nSuperS_deaths_hunt; //***b
    
        res.res_E1_hunt[tstep_0] = nE1_deaths_hunt; //6a
        res.res_SuperE1_hunt[tstep_0] = nSuperE1_deaths_hunt; //6b
    
        res.res_I_hunt[tstep_0] = nI_deaths_hunt; //7a
        res.res_SuperI_hunt[tstep_0] = nSuperI_deaths_hunt; //7b
    
        //susceptible birth/death
        res.res_nS_birth[tstep_0] = nS_births; //8a
        res.res_nSuperS_birth[tstep_0] = nSuperS_births; //8b
    
        res.res_nS_death[tstep_0] = nS_deaths; //9a
        res.res_nSuperS_death[tstep_0] = nSuperS_deaths; //9b
    
        //E1 birth/death
        res.res_nE1_death[tstep_0] = nE1_deaths; //10a
        res.res_nSuperE1_death[tstep_0] = nSuperE1_deaths; //10b
    
        //I birth/death
        res.res_nI_death[tstep_0] = nI_deaths; //11a
        res.res_nSuperI_death[tstep_0] = nSuperI_deaths; //11b
    
        // number of transitions
        res.res_trans_to_E1[tstep_0] = n_trans_to_E1; //12a
        res.res_trans_to_SuperE1[tstep_0] = n_trans_to_SuperE1; //12b
        res.res_trans_to_I[tstep_0] = n_trans_to_I; //13a
        res.res_trans_to_SuperI[tstep_0] = n_trans_to_SuperI; //13b
    
        res.res_month[tstep_0] = p.start_month; //14
        res.res_quarter[tstep_0] = p.month_to_quarter.at(p.start_month) + 1; //15
    
        res.res_Total_hunt[tstep_0] = 0;
        res.res_Total_inf[tstep_0] = nE1 + nSuperE1 + nI + nSuperI;
    
        res.res_inf_hunt[tstep_0] = 0; //16
    
        res.res_FO[tstep_0] = 0;
        res.res_FO_time[tstep_0] = 0;
        res.res_Hunt_prev[tstep_0] = 0;
    
    if(p.verbose > 2){std::cout << "task complete" << std::endl;}
    
    if(p.verbose > 0){std::cout << "beginning loop over all timesteps" << std::endl;}
    
    //start loop over all timesteps
    for(int tstep = 1; tstep <= n_tsteps; ++tstep)
    {
        if(p.verbose > 3){std::cout << "start timestep " << tstep << std::endl;}
        
        //verbose statement templates
        //if(p.verbose > 1){std::cout << "task: " << std::endl;}
        //if(p.verbose > 2){std::cout << "task complete" << std::endl;}
        
        //assign month and quarter
        if(p.verbose > 3){std::cout << "task: assign month and quarter" << std::endl;}
        
        int current_month = (tstep + p.start_month - 2) % 12 + 1;
        int current_quarter_idx = p.month_to_quarter.at(current_month);
        int current_quarter = p.month_to_quarter.at(current_month) + 1;
        double d_n = double(N);
        
        if(p.verbose > 3){std::cout << "task complete" << std::endl;}
        
        
        //initialize all temporary storage variables
        if(p.verbose > 4){std::cout << "task: Reset temporary storage variables" << std::endl;}
            
        // Birth variables
        nSuperS_births = 0;
        nS_births = 0;
        
        // S deaths variables
        nS_deaths = 0;
        nSuperS_deaths = 0;
        nS_deaths_hunt = 0;
        nSuperS_deaths_hunt = 0;
        
        // E1 death variables
        nE1_deaths = 0;
        nSuperE1_deaths = 0;
        nE1_deaths_hunt = 0;
        nSuperE1_deaths_hunt = 0;
        
        // I death variables
        nI_deaths = 0;
        nSuperI_deaths = 0;
        nI_deaths_hunt = 0;
        nSuperI_deaths_hunt = 0;
        
        // transition/transmission variables
        n_trans_to_E1 = 0;
        n_trans_to_SuperE1 = 0;
        n_trans_to_I = 0;
        n_trans_to_SuperI = 0;
        
        
        if(p.verbose > 4){std::cout << "task complete" << std::endl;}
        
        /*///////////////////////////////////////////////////////////////////////////////*/
        /*///////////////////////////////////////////////////////////////////////////////*/
        /*///////////////////////////////////////////////////////////////////////////////*/
        
        
        //compute event rates
        if(p.verbose > 4){std::cout << "task: Assign event rates" << std::endl;}
        
        //numerically compute mean monthly birth rate
        //midpoint approximation, with interval length of 1 and 'daily' step sizes
        int it = p.integrate_type;
        int n = 90; //integer days per month
        double step = 1.0 / double(n); //step interval length
        double integral_value = 0.0;  // signed area
        
        //conditional to determine how integration is performed. default is left-hand approximation
        if(it==1){
            for (int i = 0; i < n; i ++) {
                integral_value += step*p.max_birth_rate*exp( -p.synchrony*pow( cos( (M_PI/12)*( ((i+0.5)*step) + (current_month-1) - p.phase_shift) ), 2.0 ) );//midpoint approximation
            }
        }
        else if(it==3){
            for (int i = 0; i < n; i ++) {
                integral_value += step*( (p.max_birth_rate*exp( -p.synchrony*pow( cos( (M_PI/12) * (i*step + (current_month-1) - p.phase_shift) ) , 2.0 ))
                                          +
                                          p.max_birth_rate*exp( -p.synchrony*pow( cos( (M_PI/12) * ((i+1)*step + (current_month-1) - p.phase_shift) ) , 2.0 )) )
                                        / 2.0 );//trapezoid approximation
            }
        }
        else{
            for (int i = 0; i < n; i ++) {
                integral_value += step*p.max_birth_rate*exp( -p.synchrony*pow( cos( (M_PI/12) * (i*step + (current_month-1) - p.phase_shift) ) , 2.0 ));//lefthand approximation
            }
        }
        
        
        double BR = (N_d)*integral_value;

        
        
        double hunting_rate = 0.0;
        if(current_quarter == 4){hunting_rate = p.hunting_rate;}
        
        //monthly lambda
        double lambda = lambda_month[current_month-1]; //equivalent to scale_frac in CT model
        double lambda_frac = 1/lambda;
        //natural mortality
        double nat_mortality_rate = (N_d)*(p.mortality_rate + p.competition_scalar * pow((N_d / double(p.carrying_capacity)), p.competition_power));
        
        //transmission, separated by mechanism
        
        //deer-deer transmission (density dependence)
        double interspecific_trans_rate = p.deer_deer_contact_rate * p.transmission_rate * (double(nS+nSuperS)) * (double(nI+nSuperI) / (p.county_habitat_area));
    
        //transmission from infected farm contact
        double farm_trans_rate =  farm_contact_rate_v[current_quarter_idx] * p.transmission_rate * (nS+(p.SS_farm_contact_factor*nSuperS));
    
        //transition
        double transition_rate = (nE1 + nSuperE1) / draw_gamma_gsl(p.sigma1_shape, p.sigma1_scale);
        if(p.sigma1_shape == 0 || p.sigma1_scale == 0){transition_rate = 0;}
        //determine total event rate
        double tot_rate = BR + hunting_rate + nat_mortality_rate + interspecific_trans_rate + farm_trans_rate + transition_rate;
        //ordered birth, natural mortality, hunt mortality, transmission, transition
        
        if(verbose > 4){
            std::cout << "printing rates: " << std::endl;
            std::cout << "hunt rate: " << hunting_rate << std::endl;
            std::cout << "mort. rate: " << nat_mortality_rate << std::endl;
            std::cout << "Birth rate: " << BR << std::endl;
            std::cout << "transmission rate: " << interspecific_trans_rate << std::endl;
            std::cout << "spill. transmission rate: " << farm_trans_rate << std::endl;
            std::cout << "transition rate: " << transition_rate << std::endl;
        }
        
        if(p.verbose > 3){std::cout << "task complete" << std::endl;}
        
        
        
        
        
        /*//////////////////////////////////////
           Update Classes By Event Occurences
        //////////////////////////////////////*/
        
        int n_events = draw_poisson_gsl(lambda);
    
        if( (N == 0 && NSuper == 0)  ){
            n_events = 0;
        }
        
        if(n_events >= 1){
            for(int i = 1; i <= n_events; ++i){
                double U1 = draw_uniform_gsl(0.0, 1.0);
                double U2;
                //double event_cut = 0.0;
                
                if(U1 < (BR)*lambda_frac){
                    //event_cut = (BR)*lambda_frac;
                    if(verbose > 6){std::cout << "computing birth..." << std::endl;}
                    U2 = draw_uniform_gsl(0.0, 1.0);
                    if( U2 < p.prop_newborn_superSpreader ){
                        nSuperS_births += 1; nSuperS += 1;}//S_SS birth
                    else if( U2 < 1 ){
                        nS_births += 1; nS += 1; }//S birth
                    else{}
                }
                
                /*///////
                  death
                ///////*/
                
                //nat mortality
                else if(U1 < (BR + nat_mortality_rate)*lambda_frac){
                    //event_cut = (BR + nat_mortality_rate)*lambda_frac;
                    //determine class for mortality
                    if(verbose > 6){std::cout << "computing nat mortality..." << std::endl;}
                    U2 = draw_uniform_gsl(0.0, 1.0);
                    if( U2 < (double(nS) / N_d) ){
                        nS_deaths += 1; nS -= 1; }//S
                    else if( U2 < (double(nS + nE1) / N_d) ){
                        nE1_deaths += 1; nE1 -= 1; }//E1
                    else if( U2 < (double(nS + nE1 + nI) / N_d) ){
                        nI_deaths += 1; nI -= 1; }//I
                    else if( U2 < (double(nS + nE1 + nI + nSuperS) / N_d) ){
                        nSuperS_deaths += 1; nSuperS -= 1; }//S_SS
                    else if( U2 < (double(nS + nE1 + nI + nSuperS + nSuperE1) / N_d) ){
                        nSuperE1_deaths += 1; nSuperE1 -= 1; }//E1_SS
                    else if( U2 < (double(nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI) / N_d) ){
                        nSuperI_deaths += 1; nSuperI -= 1; }//I_SS
                    else{}
                }
                
                //hunt mortality
                else if(U1 < (BR + nat_mortality_rate + hunting_rate)*lambda_frac){
                    //event_cut = (BR + nat_mortality_rate + hunting_rate)*lambda_frac;
                    //determine hunted class
                    if(verbose > 6){std::cout << "computing hunt mortality..." << std::endl;}
                    U2 = draw_uniform_gsl(0.0, 1.0);
                    if( U2 < (double(nS) / N_d) ){
                        nS_deaths_hunt += 1; nS -= 1; }//S % SMS - changed nS_deaths to nS_deaths_hunt 8/15/2025
                    else if( U2 < (double(nS + nE1) / N_d) ){
                        nE1_deaths_hunt += 1; nE1 -= 1; }//E1
                    else if( U2 < (double(nS + nE1 + nI) / N_d) ){
                        nI_deaths_hunt += 1; nI -= 1; }//I
                    else if( U2 < (double(nS + nE1 + nI + nSuperS) / N_d) ){
                        nSuperS_deaths_hunt += 1; nSuperS -= 1;}//S_SS
                    else if( U2 < (double(nS + nE1 + nI + nSuperS + nSuperE1) / N_d) ){
                        nSuperE1_deaths_hunt += 1; nSuperE1 -= 1;}//E1_SS
                    else if( U2 < (double(nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI) / N_d) ){
                        nSuperI_deaths_hunt += 1; nSuperI -= 1;}//I_SS
                    else{}
                }
                
                /*/////////////////
                    Transmission
                /////////////////*/
                
                else if(U1 < (BR  + nat_mortality_rate + hunting_rate + interspecific_trans_rate)*lambda_frac){
                    //event_cut = (BR  + nat_mortality_rate + hunting_rate + interspecific_trans_rate)*lambda_frac;
                    //divide between normal and SS susceptible individuals
                    if(verbose > 6){std::cout << "computing i transmission..." << std::endl;}
                    U2 = draw_uniform_gsl(0.0, 1.0);
                    if( U2 < (double(nS)/ double(nS+nSuperS)) ){
                        n_trans_to_E1 += 1; nS -= 1; nE1 += 1; } //S becomes infected
                    else if( U2 < ( double(nS+nSuperS)/double(nS+nSuperS) ) ){
                        n_trans_to_SuperE1 += 1; nSuperS -= 1; nSuperE1 += 1; }//S_SS individual becomes infected
                    else{}
                        // //ensure other events are not evaluated
                }
                
                
                else if(U1 < (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate + farm_trans_rate)*lambda_frac){
                    //event_cut = (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate + farm_trans_rate)*lambda_frac;
                    //divide between normal and SS individuals
                    if(verbose > 6){std::cout << "computing s transmission..." << std::endl;}
                    U2 = draw_uniform_gsl(0.0, 1.0);
                    if( U2 < (double(nS) / double(nS+p.SS_farm_contact_factor*nSuperS)) ){
                        n_trans_to_E1 += 1; nS -= 1; nE1 += 1; } //S becomes infected
                    else if( U2 < double(nS+p.SS_farm_contact_factor*nSuperS) / double(nS+p.SS_farm_contact_factor*nSuperS) ){ n_trans_to_SuperE1 += 1; nSuperS -= 1; nSuperE1 += 1; } //S_SS individual becomes infected
                    else{}
                    // //ensure other events are not evaluated
                }
                                
                /*//////////////
                   Transition
                //////////////*/
        
                else if(U1 < (BR + hunting_rate + nat_mortality_rate + interspecific_trans_rate + farm_trans_rate + transition_rate)*lambda_frac){
                    //event_cut = (BR + hunting_rate + nat_mortality_rate + interspecific_trans_rate + farm_trans_rate + transition_rate)*lambda_frac;
                    //divide between normal or SS exposed individuals
                    if(verbose > 6){std::cout << "computing transition..." << std::endl;}
                    U2 = draw_uniform_gsl(0.0, 1.0);
                    if( U2 < (double(nE1) / double(nE1+nSuperE1)) ){
                        n_trans_to_I += 1; nE1 -= 1; nI += 1;} //E1 becomes infectious
                    else if( U2 < (double(nE1+nSuperE1) / double(nE1+nSuperE1)) ){
                        n_trans_to_SuperI += 1; nSuperE1 -= 1; nSuperI += 1;} //E1_SS individual becomes infectious
                    else{}
                }

                else{ if(verbose > 6){ std::cout << "No event..." << std::endl; } }
                
                //update N after event has been calculated
                N = nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI;
                N_d = double(nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI);
                
                
                if(verbose > 6){
                    std::cout << "Upper event cutoffs: " << std::endl;
                    std::cout << "U1: " << U1 << std::endl;
                    std::cout << "Birth rate: " << (BR)*lambda_frac << std::endl;
                    std::cout << "mort. rate: " << (BR + nat_mortality_rate)*lambda_frac << std::endl;
                    std::cout << "hunt rate: " << (BR + nat_mortality_rate + hunting_rate)*lambda_frac << std::endl;
                    std::cout << "transmission rate: " << (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate)*lambda_frac << std::endl;
                    std::cout << "spill. transmission rate: " << (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate + farm_trans_rate)*lambda_frac << std::endl;
                    std::cout << "transition rate: " << (BR + nat_mortality_rate + hunting_rate + interspecific_trans_rate + farm_trans_rate + transition_rate)*lambda_frac << std::endl;
                }
                
            }
        }
        
        //first ensure no negative changes
        if(nS_births < 0){std::cerr << "Negative births: " << nS_births <<  std::endl;} //ignore super spreader births: not yet calculated
        
        if(nSuperS_deaths < 0 or nS_deaths < 0 or nSuperS_deaths_hunt < 0 or nS_deaths_hunt < 0){std::cerr << "Negative S deaths: " << nS_deaths << "\t" << nSuperS_deaths << nS_deaths_hunt << "\t" << nSuperS_deaths_hunt << std::endl;}
        
        if(nSuperE1_deaths < 0 or nE1_deaths < 0 or nSuperE1_deaths_hunt < 0 or nE1_deaths_hunt < 0){std::cerr << "Negative E1 deaths: " << nE1_deaths << "\t" << nSuperE1_deaths << nE1_deaths_hunt << "\t" << nSuperE1_deaths_hunt << std::endl;}
        
        if(nSuperI_deaths < 0 or nI_deaths < 0 or nSuperI_deaths_hunt < 0 or nI_deaths_hunt < 0){std::cerr << "Negative I deaths: " << nI_deaths << "\t" << nSuperI_deaths << nI_deaths_hunt << "\t" << nSuperI_deaths_hunt << std::endl;}
        
        if(n_trans_to_E1 < 0 or n_trans_to_SuperE1 < 0 or n_trans_to_I < 0 or n_trans_to_SuperI < 0){std::cerr << "Negative transitions: " << n_trans_to_E1 << "\t" << n_trans_to_SuperE1 << n_trans_to_I << "\t" << n_trans_to_SuperI << std::endl;}
        
        
        /*apply changes
        nS =  nS - n_trans_to_E1
                 - nS_deaths - nS_deaths_hunt
                 + nS_births;
        
        nE1 = nE1 - n_trans_to_I + n_trans_to_E1
                  - nE1_deaths - nE1_deaths_hunt;
        
        nI =  nI + n_trans_to_I
                 - nI_deaths - nI_deaths_hunt;
        
        nSuperS =  nSuperS - n_trans_to_SuperE1
                           - nSuperS_deaths - nSuperS_deaths_hunt
                           + nSuperS_births;
        
        nSuperE1 = nSuperE1 - n_trans_to_SuperI + n_trans_to_SuperE1
                            - nSuperE1_deaths - nSuperE1_deaths_hunt;
        
        nSuperI =  nSuperI + n_trans_to_SuperI
                           - nSuperI_deaths - nSuperI_deaths_hunt;
    */
        
        N = nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI;
        N_d = double(nS + nE1 + nI + nSuperS + nSuperE1 + nSuperI);
        NSuper = nSuperS + nSuperE1 + nSuperI;
        
        /*/////////////
         Prelim. check
        /////////////*/
        
        //FIX THIS fix this -- for testing
        if( nS < 0 or nE1 < 0  or nI < 0 or
            nSuperS < 0 or nSuperE1 < 0  or nSuperI < 0 )
        {
            std::cerr << "Negative n animals: " << nS << "\t" << nE1 << "\t"  << nI << "\t" << nSuperS << "\t" << nSuperE1 << "\t" << nSuperI << std::endl;
        }
        
        
        
        int tot_hunt = nS_deaths_hunt + nE1_deaths_hunt + nI_deaths_hunt + nSuperS_deaths_hunt + nSuperE1_deaths_hunt + nSuperI_deaths_hunt;
        int inf_hunt = nE1_deaths_hunt + nI_deaths_hunt + nSuperE1_deaths_hunt + nSuperI_deaths_hunt;
        int n_inf = nE1 + nSuperE1 + nI + nSuperI;
        
        if(n_inf == 0){
            if(FO==0){FO_time=tstep;}
            FO = 1;
        }
        
        double hunt_prevalence = ( double(inf_hunt)/double(tot_hunt) ) / ( double(n_inf) / (N_d) );
        if( N_d == 0 or tot_hunt == 0 or n_inf == 0){ hunt_prevalence = 0; }
        //if( n_inf == 0 and inf_hunt > 0 ){ ---- }
        
        /*/////////////////////////
                Save Results
        /////////////////////////*/
        if(p.verbose > 1){std::cout << "begin saving results" << std::endl;}
            res.res_tstep[tstep] = tstep; //1
        
            res.res_n[tstep] = N; //2a
            res.res_nSuper[tstep] = NSuper; //2b
        
            res.res_nS[tstep] = nS; //3a
            res.res_nSuperS[tstep] = nSuperS; //3b
        
            res.res_nE1[tstep] = nE1; //4a
            res.res_nSuperE1[tstep] = nSuperE1; //4b
        
            res.res_nI[tstep] = nI; //5a
            res.res_nSuperI[tstep] = nSuperI; //5b
        
            res.res_S_hunt[tstep] = nS_deaths_hunt; //***a
            res.res_SuperS_hunt[tstep] = nSuperS_deaths_hunt; //***b
        
            res.res_E1_hunt[tstep] = nE1_deaths_hunt; //6a
            res.res_SuperE1_hunt[tstep] = nSuperE1_deaths_hunt; //6b
        
            res.res_I_hunt[tstep] = nI_deaths_hunt; //7a
            res.res_SuperI_hunt[tstep] = nSuperI_deaths_hunt; //7b
        
            //susceptible birth/death
            res.res_nS_birth[tstep] = nS_births; //8a
            res.res_nSuperS_birth[tstep] = nSuperS_births; //8b
        
            res.res_nS_death[tstep] = nS_deaths; //9a
            res.res_nSuperS_death[tstep] = nSuperS_deaths; //9b
        
            //E1 birth/death
            res.res_nE1_death[tstep] = nE1_deaths; //10a
            res.res_nSuperE1_death[tstep] = nSuperE1_deaths; //10b
        
            //I birth/death
            res.res_nI_death[tstep] = nI_deaths; //11a
            res.res_nSuperI_death[tstep] = nSuperI_deaths; //11b
        
            // number of transitions
            res.res_trans_to_E1[tstep] = n_trans_to_E1; //12a
            res.res_trans_to_SuperE1[tstep] = n_trans_to_SuperE1; //12b
            res.res_trans_to_I[tstep] = n_trans_to_I; //13a
            res.res_trans_to_SuperI[tstep] = n_trans_to_SuperI; //13b
        
            res.res_month[tstep] = current_month; //14
            res.res_quarter[tstep] = current_quarter_idx + 1; //15
        
            res.res_Total_hunt[tstep] = tot_hunt;
            res.res_Total_inf[tstep] = nE1 + nSuperE1 + nI + nSuperI;
            res.res_inf_hunt[tstep] = nE1_deaths_hunt + nI_deaths_hunt + nSuperE1_deaths_hunt + nSuperI_deaths_hunt; //16
        
            res.res_FO[tstep] = FO;
            res.res_FO_time[tstep] = FO_time;
            res.res_Hunt_prev[tstep] = hunt_prevalence;
        
        if(p.verbose > 2){std::cout << "results saved" <<  std::endl;}
        
        if(p.verbose > 3){std::cout << "end timestep " << tstep << std::endl;}
    }//end loop over timesteps
    
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
    std::stringstream bs;
    
    
    //Troubleshooting parameter values from wrapper code
    if(p.verbose > 4){
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
    
    if(p.verbose == -19){
        //simulation peoperties
        bs  << "rep;" << "tstep;"
            //General population properties
            << "N;"
            //misc. properties
            << "Total Hunt;"
            << "Total Infected;"
            << "Infected Hunted;"
            << "quarter;"
            << "fadeout;"
            << "fadeout time;"
            << "Hunt Prevalence"
            << std::endl;
        *out_stream << bs.rdbuf();
        
        for(int rep_i=0; rep_i<n_reps; ++rep_i){
            
            ReplicateResult res = btb_wl_model_DTMC(p);
                
            std::stringstream rep_bs;
            
            for(int i=0; i<res.n_tsteps; ++i){
                    
                rep_bs << rep_i+1 << ";" //0*
                    << i << ";" //1
                    << res.res_n[i] << ";" //2a
                    << res.res_Total_hunt[i] << ";"
                    << res.res_Total_inf[i] << ";"
                    << res.res_inf_hunt[i] << ";"
                    << res.res_quarter[i] << ";"
                    << res.res_FO[i] << ";"
                    << res.res_FO_time[i] << ";"
                    << res.res_Hunt_prev[i]
                    << std::endl;
                        
                }
            
            *out_stream << rep_bs.rdbuf();
            
        }
    }
    else{
        //simulation peoperties
        bs  << "rep;" << "tstep;"
            //General population properties
            << "N;" << "S;" << "E1;" << "I;" << "S hunted;" << "E1 hunted;" << "I hunted;"
            << "S Births;" << "S Deaths;" << "E1 Deaths;" << "I Deaths;"
            << "Trans. to E1;" << "Trans. to I;"
            //Super spreader population properties
            << "N Super;" << "SS_S;" << "SS_E1;" << "SS_I;" << "SS_S hunted;" << "SS_E1 hunted;" << "SS_I hunted;"
            << "SS_S Births;" << "SS_S Deaths;" << "SS_E1 Deaths;" << "SS_I Deaths;"
            << "Trans. to SS_E1;" << "Trans. to SS_I;"
            //misc. properties
            << "Total Hunt;"
            << "Total Infected;"
            << "Infected Hunted;" //was birth rate
            << "quarter;"
            << "fadeout;"
            << "fadeout time;"
            << "Hunt Prevalence"
            << std::endl;
        *out_stream << bs.rdbuf();
    
        for(int rep_i=0; rep_i<n_reps; ++rep_i){
            
            ReplicateResult res = btb_wl_model_DTMC(p);
                
            std::stringstream rep_bs;
            
            //only print out last timestep for sensitivity runs
            if(p.verbose == -3){
                int t_end = res.n_tsteps - 1;
                rep_bs << rep_i+1 << ";" //0*
                << t_end << ";" //1
                
                << res.res_n[t_end] << ";" //2a
                << res.res_nS[t_end] << ";" //3a
                << res.res_nE1[t_end] << ";"//4a
                << res.res_nI[t_end] << ";" //5a
                << res.res_S_hunt[t_end] << ";" //***a
                << res.res_E1_hunt[t_end] << ";" //6a
                << res.res_I_hunt[t_end] << ";" //7a
                << res.res_nS_birth[t_end] << ";" //8a
                << res.res_nS_death[t_end] << ";" //9a
                << res.res_nE1_death[t_end] << ";" //10a
                << res.res_nI_death[t_end] << ";" //11a
                << res.res_trans_to_E1[t_end] << ";" //12a
                << res.res_trans_to_I[t_end] << ";" //13a
                
                << res.res_nSuper[t_end] << ";" //2b
                << res.res_nSuperS[t_end] << ";" //3b
                << res.res_nSuperE1[t_end] << ";"//4b
                << res.res_nSuperI[t_end] << ";" //5b
                << res.res_SuperS_hunt[t_end] << ";" //***b
                << res.res_SuperE1_hunt[t_end] << ";" //6b
                << res.res_SuperI_hunt[t_end] << ";" //7b
                << res.res_nSuperS_birth[t_end] << ";" //8b
                << res.res_nSuperS_death[t_end] << ";" //9b
                << res.res_nSuperE1_death[t_end] << ";" //10b
                << res.res_nSuperI_death[t_end] << ";" //11b
                << res.res_trans_to_SuperE1[t_end] << ";" //12b
                << res.res_trans_to_SuperI[t_end] << ";" //13b
                    
                << res.res_Total_hunt[t_end] << ";"
                << res.res_Total_inf[t_end] << ";"
                << res.res_inf_hunt[t_end] << ";"
                << res.res_quarter[t_end] << ";"
                << res.res_FO[t_end] << ";"
                << res.res_FO_time[t_end] << ";"
                << res.res_Hunt_prev[t_end]
                << std::endl;
            }
            //otherwise print entire replicate as normal
            else{
                for(int i=0; i<res.n_tsteps; ++i){
                    
                    rep_bs << rep_i+1 << ";" //0*
                    << i << ";" //1
                    
                    << res.res_n[i] << ";" //2a
                    << res.res_nS[i] << ";" //3a
                    << res.res_nE1[i] << ";"//4a
                    << res.res_nI[i] << ";" //5a
                    << res.res_S_hunt[i] << ";" //***a
                    << res.res_E1_hunt[i] << ";" //6a
                    << res.res_I_hunt[i] << ";" //7a
                    << res.res_nS_birth[i] << ";" //8a
                    << res.res_nS_death[i] << ";" //9a
                    << res.res_nE1_death[i] << ";" //10a
                    << res.res_nI_death[i] << ";" //11a
                    << res.res_trans_to_E1[i] << ";" //12a
                    << res.res_trans_to_I[i] << ";" //13a
                    
                    << res.res_nSuper[i] << ";" //2b
                    << res.res_nSuperS[i] << ";" //3b
                    << res.res_nSuperE1[i] << ";"//4b
                    << res.res_nSuperI[i] << ";" //5b
                    << res.res_SuperS_hunt[i] << ";" //***b
                    << res.res_SuperE1_hunt[i] << ";" //6b
                    << res.res_SuperI_hunt[i] << ";" //7b
                    << res.res_nSuperS_birth[i] << ";" //8b
                    << res.res_nSuperS_death[i] << ";" //9b
                    << res.res_nSuperE1_death[i] << ";" //10b
                    << res.res_nSuperI_death[i] << ";" //11b
                    << res.res_trans_to_SuperE1[i] << ";" //12b
                    << res.res_trans_to_SuperI[i] << ";" //13b
                        
                    << res.res_Total_hunt[i] << ";"
                    << res.res_Total_inf[i] << ";"
                    << res.res_inf_hunt[i] << ";"
                    << res.res_quarter[i] << ";"
                    << res.res_FO[i] << ";"
                    << res.res_FO_time[i] << ";"
                    << res.res_Hunt_prev[i]
                    << std::endl;
                        
                }
            }
            *out_stream << rep_bs.rdbuf();
            
        }
        
    }
                    
    return 0;
                    
}


