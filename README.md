# Summary

These sources discuss **setting up and simulating stable N-body systems**, specifically focusing on **Plummer spheres**. The process involves **positioning particles according to a desired density distribution** (like the Plummer sphere) and then **assigning them appropriate velocities** using a distribution function based on escape velocity to ensure stability. The text also covers using **random number generators** for particle placement and velocity assignment, as well as **testing system stability** using Lagrangian radii. Finally, the sources introduce **simulating the effects of gas expulsion** on embedded star clusters, which can lead to "infant mortality" in young clusters.

## Part 1

To begin, please review the **Setting up Plummer Sphere** Python code to familiarize yourself with the parameters and initial conditions necessary to achieve the desired distribution. Each section of the code is clearly labeled with its corresponding project section and provides a detailed description of each line of code, enhancing your understanding of the project’s objectives. 

## Part 2

In the preceding section, we will proceed to utilize the **nbody.py** Python code to execute the simulation. This code generates a dynamic simulation that allows for the visualization of the cluster’s temporal evolution. However, it does not operate autonomously. From Part 1, we will extract the positions and velocities components of the Plummer sphere distribution. Subsequently, we will modify the **nbody.py** file to enable its reading from the saved table.

## Part 3

This section explores the impact of gas expulsion on embedded star clusters, focusing on the effects of varying gas fractions and removal timescales. The simulations, using a Plummer sphere model for gas and Lagrangian radii plots, investigate cluster survival outcomes and the parameter space where clusters are destroyed or survive gas expulsion. See the comments made in **Embedded Star Cluster** to better understanding.

### Additional Comments

Within this repository, you will find a laboratory report named **Lab Report 2 N-Body Sim** that serves as the primary documentation for the project’s context and background. It is imperative to recognize that the laboratory report constitutes the main work, while the repository primarily functions as a support and readability tool. Additionally, the repository contains homework tasks (**Sims_Course_Project_2,_Nbody_Sims,_Star_Clusters_&_Infant_Mortality**) that include the original assignments, their explanations, and solutions presented in the laboratory report. While reading this file is not essential, it is optional and may provide further insights.
