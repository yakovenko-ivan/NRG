# NRG Tutorial: Running Your First Simulation

This tutorial guides you through the complete process of running a simulation with the NRG software. The workflow consists of two distinct phases, each using a different executable:
1.  **Generate a Problem Setup**: Use a `package_interface` executable (e.g., for a 1D laminar flame) to create a specific case folder.
2.  **Run the Simulation**: Use the `computing_module` solver to perform the calculation from within that generated folder.

We will use the **1D Laminar Burning Velocity** test case as our example.

---

## Prerequisites

*   You have successfully built the NRG software by following the [Installation Guide](INSTALLATION.md).
*   You have built both a `package_interface` target (e.g., `package_interface_1D_laminar_velocity`) and the `computing_module` target.
*   You know the location of your compiled binaries (typically in `build/bin/` or `build/bin/Release/`).

---

## Overview: The Two-Step Workflow

A key concept in NRG is the separation between **problem setup** and **solver execution**. They are distinct steps run from different directories with different executables.

1.  **Step 1 (Setup):** You run a `package_interface` program. It **must** find a `task_setup/` folder in its current working directory. It reads this template, applies specific parameters (like chemistry model, mesh resolution), and generates a new, self-contained problem folder right where you run it.
2.  **Step 2 (Solve):** You navigate into the **newly generated problem folder**. You then run the `computing_module` solver from there. It reads the configuration files in that folder and writes results to its `data_save/` subdirectory.

**Visual Workflow:**

```text
NRG Source Tree
|
| (Copy template)
v
Binary Directory (e.g., build/bin/)
|--> Contains: package_interface_*, computing_module
|
| (Execute package_interface_*)
v
Generated Problem Folder (e.g., 1D_LBV_test/.../)
|--> Contains: data_save/, data_output/, task_setup/
|
| (Execute computing_module from inside this folder)
v
Results written to data_save/
```

---

## Part 1: Generating the Problem Setup

The `package_interface` executable looks for a `task_setup/` folder **in its immediate working directory**. You must first copy the template there.

### **Step-by-Step Instructions**

#### **For Command-Line Users (All Platforms)**

1.  **Navigate to your binary directory and copy the template.**
    Run these commands from your **NRG source directory**.
    ```bash
    # Go to where your binaries are (adjust path if using Visual Studio generator)
    cd build/bin/
    # Copy the essential template folder from the source tree
    cp -r ../../package_interface/task_setup/ .
    ```
    *Now, `build/bin/` contains both `package_interface_1D_laminar_velocity` and the `task_setup/` folder.*

2.  **Run the interface executable.**
    ```bash
    # For the 1D laminar velocity case
    ./package_interface_1D_laminar_velocity
    # On Windows Command Prompt, omit the ./
    # package_interface_1D_laminar_velocity.exe
    ```

#### **For Visual Studio IDE Users (Windows ifx)**

When running from the Visual Studio GUI, the "working directory" is the folder containing the project file (`.vfproj`).

1.  **Locate the `package_interface` project directory.** It is typically inside your `build` folder (e.g., `build/package_interface/`).
2.  **Copy the template `task_setup` folder** from `NRG/package_interface/task_setup/` into this project directory.
3.  In Visual Studio, set the `package_interface` project as the startup project, choose the **Release** configuration, and run it (**Debug** will also work but be slower).

### **What Happens & What to Expect**

*   The program will run and create a **new folder** in your current working directory (e.g., `build/bin/` or the VS project dir). Its name encodes all chosen parameters, for example:
    `1D_LBV_test/nw/cartesian/FDS/KEROMNES/9.0_pcnt/dx_1.0e-04/`
*   This new folder is your **complete problem setup**. It contains:
    *   `data_save/`: **Empty initially.** This is where results will be written.
    *   `data_output/`: Contains initial field files (e.g., `fields_init.plt`).
    *   `problem_setup.log`: A text log describing the generated case.
    *   A specific `task_setup/` folder: Configured input files (`.inf`) for this case.

**Note:** You only need to copy the master `task_setup` template **once** to a given working directory. You can run `package_interface` multiple times from there to generate different cases.

---

## Part 2: Running the Simulation

The `computing_module` solver must be run **from within the specific problem folder** generated in Part 1.

### **Step-by-Step Instructions**

1.  **Navigate into the generated problem directory.**
    ```bash
    # Change into the folder that was just created.
    # The exact name was printed to the terminal at the end of Part 1.
    cd 1D_LBV_test/nw/cartesian/FDS/KEROMNES/9.0_pcnt/dx_1.0e-04/
    ```
    Verify this folder contains `data_save/`, `data_output/`, and `problem_setup.log`.

2.  **Run the solver.** Choose the simplest option:
    *   **Option A (Copy executable here):** Copy the `computing_module` executable into this problem directory, then run it.
        ```bash
        # Copy the solver from your build directory (adjust path as needed)
        cp ../../../computing_module .
        # Run it
        ./computing_module
        ```
    *   **Option B (Run from build location):** Run it directly using a relative path.
        ```bash
        ../../../build/bin/computing_module          # If using Unix Makefiles
        ../../../build/bin/Release/computing_module.exe # If using Visual Studio generator
        ```

#### **For Visual Studio IDE Users**

1.  **Copy the generated problem folder** (e.g., `1D_LBV_test/nw/...`) into the **`computing_module` project directory**.
2.  In Visual Studio, set the `computing_module` project as the startup project.
3.  **Crucial:** In the project **Properties > Debugging**, set the **Working Directory** to the path of the copied problem folder.
4.  Run the project (preferably in **Release** configuration).

---

## Part 3: Verifying a Successful Run

Check the following to confirm your simulation ran correctly:

1.  **Console Output:** The `computing_module` will print progress information (time step, simulation time, etc.). A successful run will complete without printing fatal error messages. 
2.  **Generated Files:** After the run, check:
    *   `data_save/`: Should now contain `.plt` (Tecplot format) files with names like `000000us.plt`, `000100us.plt`, etc. These files should have a size greater than zero.
    *   `problem_setup.log`: Can be opened with a text editor to review the case description.
3.  **Visualization:** The `.plt` files can be opened with scientific visualization tools like **Tecplot** or **ParaView** to inspect the flow field, temperature, species concentrations, and other results.

---

## Complete Quick-Start Example

Here is the complete command-line workflow for a standard **Release** build with `Unix Makefiles` on Linux/macOS/WSL:

```bash
# 1. Build the software
cd /path/to/NRG
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --target computing_module package_interface

# 2. Generate a 1D laminar flame setup
cd bin
cp -r ../../package_interface/task_setup/ .
./package_interface_1D_laminar_velocity

# 3. Run the simulation
cd 1D_LBV_test/nw/cartesian/FDS/KEROMNES/9.0_pcnt/dx_1.0e-04/
cp ../../../computing_module .
./computing_module

# 4. Check results
ls -la data_save/
```
---
## Troubleshooting

This section addresses common problems encountered when setting up and running NRG simulations. If your issue isn't listed here, please search for or open a new topic on the [GitHub Issues](https://github.com/yakovenko-ivan/NRG/issues) page.

| Problem & Error Message | Likely Cause | Solution |
| :--- | :--- | :--- |
| **`package_interface` crashes immediately or prints an error about a missing `task_setup` folder.** | The executable cannot find the required template directory in its current working directory. | **Copy the master `task_setup` folder** from `NRG/package_interface/task_setup/` into the **same directory** where your `package_interface_*.exe` is located (e.g., `build/bin/`). |
| **`computing_module` crashes on startup with file errors.** | The solver is being run from the wrong directory. It must be executed from inside the specific problem folder generated by `package_interface`. | **Navigate into the generated problem folder** (e.g., `1D_LBV_test/nw/.../`) before running the solver. Ensure this folder contains `data_save/` and `problem_setup.log`. |
| **No `.plt` result files are created in the `data_save/` folder.** | The simulation terminated prematurely, possibly due to a configuration error or numerical instability. | 1. Check the terminal output for error messages.<br>2. Verify the `task_setup` template was copied correctly in the first step.<br>3. Try running a **Debug** build for more detailed logs: `cmake --build . --config Debug`. |
| **"Permission denied" error when trying to run an executable (Linux/macOS).** | The binary file does not have execute permissions. | Grant execute permissions: `chmod +x package_interface_1D_laminar_velocity computing_module`. |
| **CMake or build errors during compilation.** | Issues with compiler configuration, missing dependencies, or the build environment. | Refer to the comprehensive solutions in the [Installation Guide (INSTALLATION.md)](INSTALLATION.md#troubleshooting). Common fixes include specifying the correct compiler or generator. |
| **The `--config` flag (e.g., `--config Release`) seems to have no effect when using `Unix Makefiles`.** | This is expected behavior. The `--config` flag is only for multi-config generators like Visual Studio. | For `Unix Makefiles`, you **must** specify the build type during the *configuration* step using `-DCMAKE_BUILD_TYPE=Release` (or `Debug`) in the `cmake ..` command. |

## Next Steps

Congratulations on running your first NRG simulation! Now that you understand the core workflow, you can explore the software's full capabilities. Here are several logical paths to continue your journey.

### 1. Modify and Experiment with a Case
The most direct next step is to change the parameters of your simulation.
*   **Location:** Navigate into the `task_setup/` folder within your generated problem directory (e.g., `1D_LBV_test/nw/.../task_setup/`).
*   **Action:** Open and edit the `.inf` configuration files with a text editor. You can modify physical properties, boundary conditions, numerical schemes, or output intervals.
*   **Test:** After saving your changes, run `computing_module` again from the same problem folder to see their effect.

### 2. Run a Different Type of Simulation
NRG can model various reactive flow scenarios. To run a different case (e.g., a droplet evaporation in a high-temperature environment):
1.  **Build a new interface:** Reconfigure and build a different `package_interface` target.
    ```bash
    cd build
    cmake .. -DPACKAGE_INTERFACE_SOURCE=src/tests/classic_tests/3D_droplet_evaporation.f90
    cmake --build . --target package_interface
    ```
2.  **Follow the workflow:** Repeat the two-step tutorial steps from the beginning, using the new `package_interface_3D_droplet_evaporation` executable.

### 3. Visualize and Analyze Results
The primary output is in the `data_save/` folder in Tecplot (`.plt`) format.
*   **Recommended Tools:** Use **[Tecplot](https://www.tecplot.com/)** or **[ParaView](https://www.paraview.org/)** (open-source) to open the `.plt` files.
*   **What to Explore:** Plot contours and graphs of fields like temperature (`T`), pressure (`p`), velocity (`v`), and species concentrations (`Y`) to analyze the combustion process and verify your setup.

### 4. Dive into the Theory and Models
To understand the scientific and numerical foundations of your simulations:
*   **Review the source code** and module descriptions in the `package_library/` directory.
*   **Consult academic references** for the implemented models (e.g., specific combustion mechanisms, the CABARET scheme).
*   **Future Resource:** A comprehensive **Theory Guide** will be added to the project documentation, detailing governing equations, chemical kinetics, turbulence models, and numerical methods.

### 5. Contribute to the NRG Project
Interested in developing new features or improving the code?
*   **Explore the Structure:** Familiarize yourself with the key directories: `computing_module/` (core solvers), `package_library/` (shared physics modules), and `package_interface/` (case setup).
*   **Check Guidelines:** Look for a `CONTRIBUTING.md` file in the repository for coding standards and pull request procedures.
*   **Future Resource:** A detailed **Developer Guide** will be created to explain the architecture, data structures, and how to extend the software.

### 6. Explore Advanced Capabilities
As you become more proficient, consider exploring:
*   **Parallel Computation:** Running larger cases using MPI for distributed-memory parallel processing.
*   **Custom Models:** Implementing user-defined source terms or new physical models within the modular framework.

We hope this tutorial provided a solid foundation. Your feedback and contributions are welcome as the NRG project and its documentation continue to evolve.

