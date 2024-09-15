# Overview and Purposes
This project was originally not developed by me and was intended for researching the effects of various cardiac pumps on heart function. It was largely built around the QT framework.

The project performs calculations of a one-dimensional model of the systemic circulation (the large circle of blood circulation) with the left side of heart and corresponding boundary conditions. The structure of the vascular graph, its parameters, and boundary conditions are described in [the article](https://pubmed.ncbi.nlm.nih.gov/26100764/) (refer to Figure 10).

The code, in places, is written in a confusing way. In my defense, I inherited it in this form.

Currently, the project does not require QT, and no special setup is needed. It is assembled using CMake and executed from the console. I have run it on Linux, and it should work without issues on Mac. However, there may be difficulties with Windows, but any IDE should be able to compile the project with some adjustments.

If Windows proves problematic, it is recommended to use [virtual machine](https://www.virtualbox.org/) with a lightweight Linux distribution.

The project's objective has also changed. Now, this is the code designed for parameter tuning for a specific _"virtual"_ patient. Practically, we input age, systolic and diastolic pressure, stroke volume, and heart rate, and the program uses the _Unscented Kalman Filter_ algorithm to gradually select the necessary parameters to achieve the desired outcomes.

# Compilation Instructions

After downloading, create an empty _/build_ directory at the root, which serves a technical purpose as it will store build files.

The simplest method - is to execute the `./build.sh` bash script from the terminal. Alternatively, you can enter the following command in the console:

```cmake -DCMAKE_BUILD_TYPE=Release CMakeLists.txt && cmake --build . --parallel 8```

After this, an executable file _bloodflow_ should appear in the directory.

# Execution Instructions

## Initial Execution Without Parameter Tuning

Please note that this is the console application. This makes it convenient to run on a cluster._(where it is commonly used)_
A curious reader, slightly familiar with C++, can easily find all possible execution options in the ```main.cpp``` file, located in the ```std::string mode```.
To execute, input the following in the terminal: ```./bloodflow -option_name option_param_1 option_param_2```. The ```-``` is followed by the option name, and then the required parameters are listed. There can be an arbitrary number of parameters.

Initially use the option: ```./bloodflow -detailed age HR data/back/base.csv```, where age - is the age value between 25 and 75 inclusive with a step of 10 years, and HR - is the heart rate. The next parameter is the path to the standard set of parameters on which the model was previously configured.

The program will output CSV files with data in the _"out"_ folder, where the filename corresponds to the vessel name or Windkessel surface. For each vessel, the following information is output:
```
virtual void print_info(std::ostream &os) final {
        if (T >= 9.0 * this -> Period && T <= 10.0 * this -> Period) {
            os << T << "," << this -> get_center_pressure() * bar_to_mmhg <<
                "," << this -> get_center_speed() <<
                "," << this -> get_center_area() <<
                "," << this -> get_center_flow() <<
                "," << this -> get_mean_center_flow() << std::endl;
        }
    }
```
That is, for one cycle from the ninth to the tenth, the pressure, flow velocity, lumen area, flow, and cycle-averaged flow will be recorded, separated by commas.

A separate file is created for the heart, where data is output according to the following rule:
```
void print_info(std::ostream &os) {
        if ((T > 9.0 * Period) && (T < 10.0 * Period)) {
            os << T << "," << y0(q_av) << "," << y0(p_vent) << "," << y0(v_v) << "," << y0(ao_valve) << ","
               << y0(q_mi) << "," << y0(p_atri) << "," << y0(v_a) << "," << y0(mi_valve) << std::endl;
        }
    }
```
That is, for one cycle from the ninth to the tenth, the flow through the aortic valve, pressure in the left ventricle, volume of the left ventricle, angle of the aortic valve opening, flow through the mitral valve, pressure in the left atrium, volume of the left atrium, and angle of the mitral valve opening will be recorded.

Using this data, you can use _gnuplot_ or _python_ to create various graphs. Information with the Windkessel surfaces can be found similarly, but it is not of current interest.

## Execution With Parameter Tuning
To do this, use the following option: ```./bloodflow -manual age HR SBP DBP SV HR```. Specify the parameters: age (from 25 to 75 with a step of 10), heart rate, systolic pressure, diastolic pressure, stroke volume.

Next, the parameter tuning process will start, which can be lengthy. The process is described in the ```ukf.cpp```, where the _Unscented Kalman Filter_ is implemented. Changes and updates can be observed in [the draft article](https://cloud.mail.ru/public/pxdG/eJRhL1yed).

After run, a log file will be created in the _data/log_ folder with information on the parameter tuning process.

## Other Options
There are more options, but they essentially act as wrappers over the two previous ones. The difference is that manually entering target pressures and stroke volumes each time isn't desirable. Instead, a configuration file listing everything can be created, sent to the cluster and results collected.
There is also an option related to Bayesian optimization, which is more complex as it involves a separate Python script. This will be considered later.

# Summary of Files
- _main.cpp_ - Entry point of the application.
- _ukf.h ukf.cpp_ - Class responsible for parameter tuning using the Kalman filter.
- _task.h task.cpp_ - Important class defining the computational task. It determines heart and vessel parameters from configuration JSON files.
- _matrix_utils.h matrix_utils.cpp_ - Auxiliary class for converting CSV files to matrices and in reverse way.
- _bayes_runner.h bayes_runner.cpp_ - Class for parameter tuning using Bayesian optimization.
- _vertex.h vertex.cpp_ - Class describing the "vertex" of the vascular graph, or rather some boundary conditions. Includes vascular bifurcation points and connection to "simplified heart", where instead of a separate heart model, the function Q(t) is used.
- _edge.h edge.cpp_ - Class describing the edge of the vascular graph. Includes its parameters, numerical scheme, result recording.
- _graph.h graph.cpp_ - Class connecting vessels with their vertices. Not explicitly used.
- _calculator.h calculator.cpp_ - Auxiliary class for calculating average values during one cardiac cycle.
- _rescaler.h_ - Class responsible for setting a number of parameters tuned by the Kalman filter. For example, it selects the total capacity and resistance of Windkessel elements. These need recalculation for each element and recording in the configuration. Only then can the task be correctly executed. (Yes, this is an odd implementation detail, but it's simpler to change everything in the configuration, and _Task_ class will read it.)
- _rcwindkessel.h_ - Class for describing Windkessel elements.
- _heart_part.h heart_part.cpp_ - Class responsible for describing and calculating the heart model.
- _heart_advanced_valves.h_ - Class responsible for the heart model with valves.
- _heart_reg.h_ - Class responsible for the heart model with valvular regurgitation pathology.
- _detailed_run.h_ - Class for running calculations with specified parameters and saving all data to files for further analysis.
- _csv_reader.h_ - Auxiliary class for working with CSV files.

**That's all!**
For questions, contact via Telegram: @rogovrt
