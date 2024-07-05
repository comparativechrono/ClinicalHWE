# ClinicalHWE
An app for running and interpreting HWE on clinical variants
Sure! Here's a detailed README file for your GitHub repository:

# Hardy-Weinberg Equilibrium (HWE) Analysis App

This Streamlit application allows you to perform Hardy-Weinberg Equilibrium (HWE) analysis with advanced features including:
- Allele Frequency Evolution Simulation
- Statistical Power Analysis
- Advanced Statistical Tests
- Clinical Interpretation of HWE Deviations
- Genotype-Phenotype Association Analysis

## Features

### Allele Frequency Evolution Simulation
- Simulate changes in allele frequencies over generations under different scenarios such as selection, migration, and population size changes.
- Visualize the allele frequency evolution over generations.

### Statistical Power Analysis
- Calculate the statistical power of detecting deviations from HWE given sample size, effect size, and significance level.
- Visualize the power curve to understand the relationship between sample size and power.

### Advanced Statistical Tests
- Perform the HWE chi-square test to assess whether the observed genotype frequencies deviate from the expected frequencies under HWE.
- Perform Fisher's exact test for small sample sizes.
- Visualize observed and expected genotype counts.
- Provide a clear explanation of how expected counts are calculated.

### Clinical Interpretation of HWE Deviations
- Interpret deviations from HWE in a clinical context, including potential causes such as population stratification, selection, or genotyping errors.
- Calculate and provide an HWE cutoff for clinical interpretation.

### Genotype-Phenotype Association Analysis
- Analyze the association between genotypes and phenotypes.
- Perform chi-square tests for association between genotype and phenotype.
- Fit logistic regression models to predict phenotypes based on genotype data.
- Visualize the distribution of phenotypes across different genotypes.
- Calculate and display the relative risk or odds ratio for having a particular phenotype given the genotype.

## Installation

To run this app locally, follow these steps:

1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-username/hwe_analysis_app.git
   cd hwe_analysis_app
   ```

2. **Create a virtual environment**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows use `venv\Scripts\activate`
   ```

3. **Install the required packages**:
   ```bash
   pip install -r requirements.txt
   ```

4. **Run the Streamlit app**:
   ```bash
   streamlit run hwe_app.py
   ```

## Usage

### Home Page
Navigate through the different features using the sidebar on the left. The home page provides an overview of the app.

### Allele Frequency Evolution
Adjust the parameters using the sliders and input fields, then click "Run Simulation" to visualize the evolution of allele frequencies over generations.

### Statistical Power Analysis
Input the effect size, sample size, and significance level, then click "Calculate Power" to see the statistical power. A power curve will also be displayed.

### Advanced Statistical Tests
Input the observed genotype counts and view the calculated expected counts under HWE. The app will perform the HWE chi-square test and Fisher's exact test (if requested). Visualizations and tables are provided for easier interpretation.

### Clinical Interpretation
Input the observed genotype counts and the app will help interpret deviations from HWE, providing guidelines and a calculated HWE cutoff for clinical interpretation.

### Genotype-Phenotype Association
Input the observed genotype counts and phenotype data. The app will perform chi-square tests for association and logistic regression, providing visualizations and statistical results.

## Contributing

Contributions are welcome! If you have any suggestions for new features or improvements, feel free to open an issue or submit a pull request.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

Special thanks to the contributors and the open-source community for their valuable tools and libraries.

---

For any queries or support, please contact [tjh70@cam.ac.uk](mailto:tjh70@cam.ac.uk).

```
.
