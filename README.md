# Automated-Matrix-Element-Generation

## Overview

The Automated-Matrix-Element-Generation (AMEG) project provides a robust framework for generating matrix elements automatically in particle physics calculations. It generates Feynman diagrams and computes Scattering Cross-sections for a given QED three-level interaction. It leverages Python for scripting and computation, offering tools to facilitate the complex process of calculating matrix elements essential for understanding particle interactions.

## Contents

- ameg.py: The core script that implements the automated generation of matrix elements.
- utils1.py & utils2.py: Utility modules that provide support functions and classes for mathematical operations and data handling based on[Spacetime-Calc-Engine](https://github.com/adyoussef/Spacetime-Calc-Engine) and [Quantum-Particle-Simulator](https://github.com/adyoussef/Quantum-Particle-Simulator).
- example.ipynb: A Jupyter notebook that demonstrates the use of AMEG with practical examples and tutorials.
- ParticleData.xml: An XML file containing data on particles used in the matrix element calculations.
- Documentation (ameg.pdf): A comprehensive guide and documentation for the AMEG project, detailing the theoretical background, software architecture, and user instructions.

## Getting Started

The externaml package [feynmna](https://pypi.org/project/feynman/) was used the generate the feynman diagrams, which can be installed with

```
pip install feynman
```

```
git clone https://github.com/adyoussef/Automated-Matrix-Element-Generation.git
cd Automated-Matrix-Element-Generation
```

## Usage
Refer to example.ipynb and Ameg.pdf for a demonstration and for the documentation



