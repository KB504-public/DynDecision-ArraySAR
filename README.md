# DynDecision-ArraySAR

![Overview of DynDecision-ArraySAR](framework_overview.jpg)

In this work, as shown above, we propose a dynamic decision-making framework for array-SAR imaging. Unlike conventional approaches that rely on static structures, fixed parameters, and outcome-level supervision, our method introduces dynamic structures and parameters guided by state-level supervision, enabling adaptive processing across diverse measurement conditions. The imaging process is reformulated from a static pipeline into a Markov decision process, where the decision-making module selects actions according to the current state, the state-transition module updates the state, and the evaluation module provides feedback on both states and actions. Through this perspective of "from states to sequence, and decision making," the framework advances array-SAR imaging from rigid fixed flows toward adaptive and scenario-specific imaging. This shift allows the proposed method to generalize better for large-scale sensing of diverse scenarios.

---

## 📖 How to access?

Our research team maintains specific guidelines for sharing our primary code and associated datasets. For the baseline methods, code would be avaiable without requiremetns and for the primal code,
To obtain access to the code, please complete and sign the required agreement, then email the scanned document to xdma@std.uestc.edu.cn or zhanxu@std.uestc.edu.cn.




This repository contains the **MATLAB implementation** of our paper:  

**Beyond Static Imaging: A Dynamic Decision Paradigm for Robust Array-SAR in Diverse Sensing Scenarios**

Currently, we provide the **StatFilter** implementation, denoting the static matched filtering method.  
This represents the earliest paradigm in array-SAR imaging, where the reconstruction is directly obtained via a fixed analytical operator, without leveraging any prior knowledge of the scene.


---

## ⚙️ Features
- Uses output_rawdata.mat (download from cloud link https://drive.google.com/drive/folders/1iG6086eDNRvF7BJp8Ej264ZIsaQhE70H?usp=drive_link)
- Implements StatFilter (static matched filtering paradigm)
- Generates reconstruction results via fixed analytical operator
- Supports dB visualization with configurable clipping ranges

---

## ▶️ Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/KB504-public/DynDecision-ArraySAR.git
   cd DynDecision-ArraySAR

---

## 🔜 TODO
- Add decision-making module code (state–sequence–decision framework)
- Add training and evaluation scripts
- Provide sample dataset generator for testing
- Add Python implementation 

---

## 📜 License

This project is licensed under the **Apache License 2.0**.  
You may freely use, modify, and distribute this code, provided that you comply with the terms of the license.  

See the [LICENSE](LICENSE) file for full details.


---


