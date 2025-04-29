# panco

Panco (Performance Analysis with Network Calculus and Optimization) regroups a set of modules for the implementation of deterministic network calculus results to compute worst-case performance in networks. It focuses on the implementation of methods using linear programming (LP) techniques (hence the optimization term, but also some more classical results (revisited with LP).

Several directories are proposed: 
- **descriptor** describes a network, by describing a set of flows and of servers
- **fifo** is the implementation of the fifo multiplexing analysis
- **edf** focuses on the Earliest-Deadline-First scheduling policy
- **staticpriorities** on the priorities
- **tsn** on the implementation of some TSN mechanisms, mainly the AVB/CBS scheduling (with taking into account the gate openings for the ime-triggered flows.

# References
The modules are mainly basec on the following publications:

  [1] Anne Bouillard,  Marc Boyer and Euriell Le Corronc. Deterministic Network Calculus: From Theory to Practical Implementation, . ISBN: 978-1-119-56341-9, October 2018, Wiley-ISTE, 352 pages (**Descriptor**)

  [2] Anne Bouillard and Giovanni Stea. Exact Worst-case Delay for FIFO-multiplexing tandems. Sixth International Conference on Performance Evaluation Methodologies and Tools (ValueTools 2012) (**FIFO**)
  
  [3] Anne Bouillard: Trade-off between accuracy and tractability of Network Calculus in FIFO networks. Perform. Evaluation 153: 102250 (2022) (**FIFO**)
   
  [4] Anne Bouillard. Admission Shaping With Network Calculus. IEEE Netw. Lett. 6(2): 115-118 (2024) (**FIFO-AdmTFA**)
  
  [5] Anne Bouillard, Laurent Jouhet and Eric Thierry. Tight Performance Bounds in the Worst-Case Analysis of Feed-Forward Networks, 29th IEEE International Conference on Computer Communications (INFOCOM 2010), p. 1316-1324 (**Blind**)
  
  [6] Anne Bouillard, Éric Thierry and Laurent Jouhet. Service curves in Network calculus: dos and don’ts RR INRIA 7094, 2009. (**static pririties**)
  
  [7] Luxi Zhao, Paul Pop, Zhong Zheng, Hugo Daigmorte, and Marc Boyer. Latency analysis of multiple classes of AVB traffic in TSN with standard credit behavior using network calculus. IEEE Trans. Ind. Electron., 68(10):102911--10302,  2021. (**TSN**)

  [8] Anne Bouillard. Earliest-deadline-first in sink trees, WONECA 2022. [video](https://www.youtube.com/watch?v=4B0ST5TsGiI&list=PLGrWRLGd9yS_nezfKdxK1x-e3yNt1krwj&index=14]), [slides](https://drive.google.com/file/d/12Vvblys74SuEuMzeQ5wWL7cSlr3ngy0n/view) (**EDF**)

 
  


# Requirements: 
Python 3, with numpy package installed, and [lp_solve](https://sourceforge.net/projects/lpsolve/)

# Installation:
```
pip install .
````
in the directory either globally or in a virtual environment

# Documentation: 
The *work in progress* documentation can be found [there](https://www.di.ens.fr/~bouillar/Panco/html/index.html)
