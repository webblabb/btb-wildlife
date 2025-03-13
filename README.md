# bTB Wildlife Models Repository

## Introduction
This repository contains code to run bTB (bovine Tuberculosis) Wildlife Models and generate figures for the manuscript titled "Stochastic modeling of bovine tuberculosis dynamics in white-tailed deer".

### C++ Compilation
Ensure you have a C++ compiler installed. You can compile the C++ files using the provided Makefile or manually using the following commands:
```sh
g++ -o wl_model_CTMC bTB_wildlifeModel_CTMC.cpp
g++ -o wl_model_DTMC bTB_wildlifeModel_DTMC.cpp
```

## Usage

### Running the Models
1. Run the R scripts to preprocess data and set up the environment:
    ```sh
    Rscript bTB_wildlife_modelPipeline.R
    ```

2. Execute the compiled C++ models:
    ```sh
    ./wl_model_CTMC
    ./wl_model_DTMC
    ```

3. Generate figures and analyze results using the R scripts:
    ```sh
    Rscript tb_wildlife_pipeline_v1.R
    ```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
For questions or feedback, please contact:
- [Lindsay Beck-Johnson](mailto:L.Beck-Johnson@colostate.edu)
