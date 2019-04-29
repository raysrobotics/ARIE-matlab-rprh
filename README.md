# Attractive Region in Environment for round-peg round-hole insertion

<p>
    <a href="https://github.com/raysrobotics/ARIE-matlab-rprh/blob/master/LICENSE"><img alt="Software License" src="https://img.shields.io/badge/license-CASIA-blue.svg"></a>
    <a><img alt="MATLAB Version" src="https://img.shields.io/badge/matlab-v2018-yellow.svg"></a>
</p>

These MATLAB scripts are used to visualize 3D attractive region in environment for round-peg-round-hole assembly tasks.

Author: Rui Li (raysworld@outlook.com)



### [Functions]

| Filename                 | Purpose                                                      | Usage                                                        |
| ------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| rot.m                    | Perform general rotation transformation                      | T = rot(f, theta)                                            |
| trans.m                  | Perform translation transformation                           | T = trans(vector)                                            |
| ARIE_log.m               | Write the current parameters, options and flag in a .txt log file | ARIE_log(param, options, flag)                               |
| ARIE_output.m            | Write the output point data in a .txt file                   | ARIE_output(x, y, z)                                         |
| localMaximum.m           | This function returns the indexes\subscripts of local maximum in the data x. x can be a vector or a matrix of any dimension. | varargout = localMaximum(x,minDist, exculdeEqualPoints)      |
| ARIE_findLowestPoint.m   | Find the lowest point in the given region. The radius of the peg and the hole are treated as the same in this version. | [flag, d_lowest, O_pe_H] = ARIE_findLowestPoint(param, options, display) |
| ARIE_findLowestPoint_R.m | In this implementation, the radius of the peg and the hole are different. |                                                              |
| ARIE_findLowestPoint_N.m | In this implementation, the intersection points are calculated with numerical methods instead of solving quadratic equations |                                                              |



[Scripts]
---------------------------------------------------------------------
**script_ecc_peg_hole.m**

It is a script-version of the file ARIE_findLowestPoint and you can read and test this file to know how ARIE_findLowestPoint works.

You can use this script to test the function *ARIE_findLowestPoint*. Given the pose of the peg, you will get the position of the lowest point on the peg.

---------------------------------------------------------------------
**script_ecc_peg_hole_edge.m**

It is a script-version of the file ARIE_findLowestPoint and you can read and test this file to know how ARIE_findLowestPoint works.

You can use this script to test the function *ARIE_findLowestPoint*. Given the pose of the peg, you will get the position of the lowest point on the peg. The difference between this file and the above file is that in this implementation, the contact points on the side surface are also detected.



---------------------------------------------------------------------
**script_main.m**

Basic version of the script. Given the tilt angle of the peg, the x and y range of the  position of the reference point attached to the peg, the radius of the peg and the hole, a 3D constrained region (and possibly an attractive region in environment) will be extracted. 

Note for this version, you can only specify the tilt angle about x axis OR y axis. Specifying both axes at the same time will lead to inaccurate results.

---------------------------------------------------------------------
**script_main_AR1.m**

This file is identical to the basic version *script_main*. The difference is that all Nan value are reassigned as 0 for postprocessing and visualization.

---------------------------------------------------------------------
**script_main_o.m**

This file is identical to the basic version *script_main*. The difference is that the output is the position of the center point of the bottom surface instead of the lowest point.

------

**script_main_plot.m**

This file is identical to the basic version *script_main*. The difference is that this version visualizes the region with Matlab built-in function *plot3* instead of *mesh*.

---------------------------------------------------------------------
**script_main_rot.m**

This file is identical to the basic version *script_main*. The difference is that you can specify both x and y tile angles with this version.

------

**script_main_AR1.m**

This file is identical to the basic version *script_main*. The difference is when lowest_point == NaN, the value will be set to 0.

------

**script_main_numeric.m**

This file is identical to the basic version *script_main*. The difference is that the intersection points are calculated with numerical methods instead of solving quadratic equations.



## [Cite]

If you use this toolbox in your research please cite it:

```
@article{rui2017,
   author = {Li, Rui and Qiao, Hong},
   title = {Condition and Strategy Analysis for Assembly Based on Attractive Region in Environment},
   journal = {IEEE/ASME Transactions on Mechatronics},
   volume = {22},
   number = {5},
   pages = {2218-2228},
   ISSN = {1083-4435},
   DOI = {10.1109/TMECH.2017.2705180},
   year = {2017}
}
```