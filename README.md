# DynDecision-ArraySAR

![Overview of DynDecision-ArraySAR](framework_overview.jpg)

Dynamic Decision Array-SAR reformulates 3D radar imaging as a Markov decision process, 
replacing static reconstruction with a state–sequence–decision framework. 
Adaptive actions and evaluation-guided feedback enable robust, scenario-adaptive imaging 
across noise levels, measurement models, and scene distributions.

---

## 📖 Description
This repository contains the **MATLAB implementation** of our paper:  

**Beyond Static Imaging: A Dynamic Decision Paradigm for Robust Array-SAR in Diverse Sensing Scenarios**

Currently, we provide the **2-D RMA Imaging with Multi-Range Maximum Intensity Projection (MIP)** script.  
The script demonstrates how to perform motion correction, RMA focusing, and MIP-based visualization.

---

## ⚙️ Features
- Reads raw data from `output_rawdata.mat`
- Applies motion correction in the array plane
- Implements RMA focusing with frequency-domain phase compensation
- Produces **512×512 SAR images**
- Supports dB visualization with configurable clipping ranges

---

## ▶️ Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/DynDecision-ArraySAR.git
   cd DynDecision-ArraySAR
