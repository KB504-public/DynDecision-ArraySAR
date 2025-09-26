# Nearfield Sample Data

This folder contains the measurement data for the near-field sensing test. The data corresponds to measurements of a **wrench target** using a **ground-based prototype system**.

## Target Image

For a visual reference of the target, you can view the optical image of the **wrench target** in the same directory as this README file. The image is located at:  
[nearfield_wrench.jpg](nearfield_wrench.jpg)

## Data Details

- **File Name**: `output_rawdata.mat`
- **Measurement Setup**: Ground-based prototype system
- **Target**: Wrench
- **Data Description**: This dataset includes the raw measurements collected from the system, which are used to test the `stat_filter_nearfield.m` function.

## How to Use

1. **Download the Sample Data**:  
   You can download the `output_rawdata.mat` from the following Google Drive link:  
   [Download Nearfield Sample Data](<https://drive.google.com/file/d/1QYB5DZnyG1y4DZuP6iV7tkJSdx-hpy63/view?usp=sharing>)

2. **Place the Data File**:  
   Once the file has been downloaded, place the `output_rawdata.mat` in the same directory as `stat_filter_nearfield.m`.

3. **Run the `stat_filter_nearfield.m`**:  
   You can now run the `stat_filter_nearfield.m` script to test the method with this data. The data is ready for use without additional preprocessing.

