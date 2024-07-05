import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare, fisher_exact, chi2_contingency
from statsmodels.stats.power import GofChisquarePower
import requests

# Title and introduction
st.title("Hardy-Weinberg Equilibrium Analysis App for Clinical Genetics")

# Sidebar for navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Allele Frequency Evolution", "Statistical Power Analysis", "Advanced Statistical Tests", "Clinical Interpretation"])

# Home Page
if page == "Home":
    st.header("Welcome to the HWE Analysis App")
    st.write("Use the navigation on the left to access different features.")
    st.write("""
    This app allows you to perform Hardy-Weinberg Equilibrium (HWE) analysis with advanced features including:
    - Allele Frequency Evolution Simulation
    - Statistical Power Analysis
    - Advanced Statistical Tests
    - Clinical Interpretation of HWE Deviations
    """)
    
# Function to calculate expected counts under HWE
def hwe_expected(obs_hom_var, obs_het, obs_hom_ref):
    total = obs_hom_var + obs_het + obs_hom_ref
    p = (2*obs_hom_var + obs_het) / (2*total)
    q = 1 - p
    exp_hom_var = p**2 * total
    exp_het = 2 * p * q * total
    exp_hom_ref = q**2 * total
    return exp_hom_var, exp_het, exp_hom_ref

# Function for allele frequency simulation
def simulate_allele_freq(p0, selection_coef, migration_rate, population_size, generations):
    freqs = [p0]
    p = p0
    for _ in range(generations):
        q = 1 - p
        p_next = p * (1 - selection_coef) + migration_rate * q
        p_next /= (p_next + q * (1 - selection_coef))
        freqs.append(p_next)
        p = p_next
    return freqs

# Function for power calculation
def calculate_power(effect_size, nobs, alpha):
    power_analysis = GofChisquarePower()
    power = power_analysis.solve_power(effect_size=effect_size, nobs=nobs, alpha=alpha)
    return power

# Function for exact test
def exact_test(obs):
    _, p_value = fisher_exact(obs)
    return p_value

# Function for chi-square test with multiple alleles
def chi_square_test(obs):
    chi2, p_value, _, _ = chi2_contingency(obs)
    return chi2, p_value

# Function to fetch gene and variant information
def fetch_variant_info(variant):
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{variant}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        return None

if page == "Allele Frequency Evolution":
    st.header("Allele Frequency Evolution Simulation")
    
    # User inputs
    p0 = st.slider("Initial Allele Frequency (p0)", 0.0, 1.0, 0.5)
    selection_coef = st.slider("Selection Coefficient (s)", 0.0, 1.0, 0.1)
    migration_rate = st.slider("Migration Rate (m)", 0.0, 1.0, 0.01)
    population_size = st.number_input("Population Size (N)", 1, 100000, 1000)
    generations = st.number_input("Number of Generations", 1, 1000, 100)
    
    # Run simulation
    if st.button("Run Simulation"):
        freqs = simulate_allele_freq(p0, selection_coef, migration_rate, population_size, generations)
        
        # Plot results
        plt.figure(figsize=(10, 6))
        plt.plot(freqs, label='Allele Frequency')
        plt.xlabel("Generations")
        plt.ylabel("Allele Frequency")
        plt.title("Allele Frequency Evolution Over Generations")
        plt.legend()
        st.pyplot(plt)

if page == "Statistical Power Analysis":
    st.header("Statistical Power Analysis")
    
    # User inputs
    effect_size = st.number_input("Effect Size (w)", 0.0, 2.0, 0.5)
    nobs = st.number_input("Sample Size (N)", 1, 100000, 1000)
    alpha = st.slider("Significance Level (alpha)", 0.0, 1.0, 0.05)
    
    # Calculate power
    if st.button("Calculate Power"):
        power = calculate_power(effect_size, nobs, alpha)
        st.write(f"Statistical Power: {power:.4f}")
        
        # Provide recommendations
        if power < 0.8:
            st.write("The statistical power is less than 0.8. Consider increasing the sample size for more reliable results.")
        else:
            st.write("The statistical power is sufficient for reliable results.")

        # Visualize power analysis
        sample_sizes = np.arange(10, nobs * 2, 10)
        powers = [calculate_power(effect_size, size, alpha) for size in sample_sizes]

        plt.figure(figsize=(10, 6))
        plt.plot(sample_sizes, powers, label='Power Curve')
        plt.axhline(0.8, color='r', linestyle='--', label='0.8 Power Threshold')
        plt.xlabel("Sample Size")
        plt.ylabel("Statistical Power")
        plt.title("Power Analysis Curve")
        plt.legend()
        st.pyplot(plt)

if page == "Advanced Statistical Tests":
    st.header("Advanced Statistical Tests")
    
    # User inputs
    st.write("Enter the observed genotype counts:")
    homozygous_variant = st.number_input("Homozygous Variant", 1, 100000, 1)
    heterozygous = st.number_input("Heterozygous", 1, 100000, 1)
    homozygous_reference = st.number_input("Homozygous Reference", 1, 100000, 1)
    
    obs = [homozygous_variant, heterozygous, homozygous_reference]
    exp = hwe_expected(*obs)
    
    try:
        # Perform HWE chi-square test
        chi2, p_value = chisquare(obs, exp)
        st.write(f"HWE Chi-Square Test: chi2 = {chi2:.4f}, p-value = {p_value:.4f}")
    except Exception as e:
        st.write(f"An error occurred during HWE chi-square test: {e}")
    
    try:
        # Perform exact test
        if st.button("Perform Exact Test"):
            p_value_exact = exact_test([[homozygous_variant, heterozygous], [heterozygous, homozygous_reference]])
            st.write(f"Exact Test P-Value: {p_value_exact:.4f}")
    except Exception as e:
        st.write(f"An error occurred during exact test: {e}")
    
    try:
        # Perform chi-square test for multiple alleles
        if st.button("Perform Chi-Square Test"):
            chi2_mult, p_value_mult = chi_square_test([[homozygous_variant, heterozygous], [heterozygous, homozygous_reference]])
            st.write(f"Chi-Square Statistic: {chi2_mult:.4f}")
            st.write(f"Chi-Square Test P-Value: {p_value_mult:.4f}")
    except Exception as e:
        st.write(f"An error occurred during chi-square test for multiple alleles: {e}")

if page == "Clinical Interpretation":
    st.header("Clinical Interpretation of HWE Deviations")
    
    st.write("""
    In clinical genetics, deviations from Hardy-Weinberg Equilibrium (HWE) in a cohort of patients with a specific condition can indicate several things:
    - **Population Stratification**: Differences in allele frequencies due to population substructure.
    - **Selection**: The condition may be associated with certain genotypes, leading to an overrepresentation of these genotypes in the cohort.
    - **Genotyping Errors**: Technical errors in genotyping can lead to deviations.
    """)
    
    st.write("To help interpret deviations from HWE, use the tools below to analyze your data:")
    
    # User inputs for genotype counts
    st.write("Enter the observed genotype counts:")
    homozygous_variant = st.number_input("Homozygous Variant", 0, 100000, 0, key="ci_hom_var")
    heterozygous = st.number_input("Heterozygous", 0, 100000, 0, key="ci_het")
    homozygous_reference = st.number_input("Homozygous Reference", 0, 100000, 0, key="ci_hom_ref")
    
    obs = [homozygous_variant, heterozygous, homozygous_reference]
    exp = hwe_expected(*obs)
    
    # Perform HWE chi-square test
    chi2, p_value = chisquare(obs, exp)
    st.write(f"HWE Chi-Square Test: chi2 = {chi2:.4f}, p-value = {p_value:.4f}")
    
    # Interpretation guidelines
    if p_value < 0.05:
        st.write("The p-value is less than 0.05, indicating a significant deviation from HWE.")
        st.write("This could suggest population stratification, selection for certain genotypes, or genotyping errors.")
        if homozygous_variant > exp[0]:
            st.write("The observed number of homozygous variants is higher than expected, which might suggest selection for this genotype in the cohort.")
        else:
            st.write("The observed number of homozygous variants is not higher than expected based on HWE.")
    else:
        st.write("The p-value is greater than 0.05, indicating no significant deviation from HWE.")
        st.write("The observed genotype frequencies are consistent with those expected under random mating conditions.")
    
    # Provide HWE cutoff
    st.write("### HWE Cutoff for Clinical Interpretation")
    cutoff_hom_var = exp[0] * 1.5  # Arbitrary cutoff value; adjust as needed
    st.write(f"If the observed number of homozygous variants is greater than {cutoff_hom_var:.2f}, it may suggest selection for this genotype in the cohort.")
    
    # Variant information
    variant = st.text_input("Enter Variant ID (e.g., rs12345) for Additional Information:")
    if variant:
        info = fetch_variant_info(variant)
        if info:
            st.write(f"Variant Information for {variant}:")
            st.write(info)
        else:
            st.write(f"Information for {variant} not found.")
