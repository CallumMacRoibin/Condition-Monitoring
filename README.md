# Advanced Condition Monitoring System for Rotary Machines

## Project Overview
This research project focuses on the development of an advanced condition monitoring (CM) system tailored for rotary machines, utilizing MATLAB for data analysis and machine learning-based fault detection.

Key highlights of the project include:
- Collection of a robust dataset of vibrational acceleration measurements from an induction motor test rig simulating various fault conditions (imbalance, misalignment, debris contamination, and lubrication issues).
- Development of a Convolutional Neural Network (CNN) achieving a fault classification accuracy of **97.42%**.
- Integration of advanced signal processing techniques (e.g., spectral kurtosis, band-pass filtering, and continuous wavelet transforms) to extract fault features and provide actionable insights to machine operators.
- Scalability and adaptability, making the system suitable for diverse condition monitoring applications.

---

## Project Structure
### Code Flow
The project follows a structured flow for processing data and training the model:
1. **Data Preprocessing**: Raw vibration signals are processed using spectral kurtosis and band-pass filtering.
2. **Feature Extraction**: Continuous wavelet transforms (CWT) generate time-frequency scalograms.
3. **Model Training**: CNN is trained on processed scalograms for fault classification.
4. **Result Visualization**: Annotated visuals are generated to highlight key fault indicators for actionable insights.

### Code Features
- MATLAB scripts for raw data preprocessing, feature extraction, and model evaluation.
- CNN implementation for fault classification using scalograms.
- Visualization of health indicators for predictive maintenance.

### Hardware
- **Induction Motor Test Rig**: Used for simulating fault conditions and collecting vibration data.
- **Vibration Sensor**: Captures acceleration measurements for analysis.

### Requirements
- MATLAB R2023a (or later)
- Deep Learning Toolbox
- Signal Processing Toolbox
- Hardware setup for data collection (induction motor and vibration sensor)

---

## Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/CallumMacRoibin/Condition-Monitoring.git
   cd Condition-Monitoring
2. Load your dataset into MATLAB using the provided scripts.
3. Run preprocessing.m to preprocess raw vibration data.
4. Execute train_CNN.m to train the fault classification model.
5. Visualize results using plot_results.m.

 ## Results & Discussion
The CNN demonstrated a high fault classification accuracy of **97.42%** on the test set. Key findings include:
- **Fault-specific insights**: The system accurately distinguishes between imbalance, misalignment, debris contamination, and lubrication issues.
- **Scalability**: The model generalizes well to new data, enabling use in different industrial setups.

---

## Future Work
- Expand the dataset to include additional fault types and operating conditions.
- Optimize the CNN architecture for faster inference and real-time fault detection.
- Develop a user-friendly interface for non-technical operators.

---

## Conclusion
This project demonstrates the effectiveness of an advanced CM system for rotary machines, combining MATLAB-based signal processing with machine learning. The solution offers a scalable and adaptable framework for predictive maintenance across various industrial applications.

