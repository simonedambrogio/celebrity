# Celebrity Endorsement and Visual Attention Study

This repository contains the data and analysis code for the research paper:

> D'Ambrogio, Werksman, Platt, & Johnson (2023). How celebrity status and gaze direction in ads drive visual attention to shape consumer decisions. *Psychology & Marketing*, 40(4), 741-755. [DOI: 10.1002/mar.21772](https://onlinelibrary.wiley.com/doi/full/10.1002/mar.21772)

## Repository Structure

```
.
├── Data/
│   ├── data_ddm.csv         # Drift-diffusion model and logistic regression data
│   ├── data_fixations.csv   # Eye fixation data from advertisement phase
│   ├── data_pupil.csv       # Processed pupil diameter data
│   └── README.txt           # Detailed data description
├── Stan/
│   ├── hDDM1.stan          # Hierarchical Drift Diffusion Model 1
│   ├── hDDM2.stan          # Hierarchical Drift Diffusion Model 2
│   ├── hDDM3.stan          # Hierarchical Drift Diffusion Model 3
│   ├── haDDM2.stan         # Hierarchical Attentional Drift Diffusion Model 2
│   └── haDDM3.stan         # Hierarchical Attentional Drift Diffusion Model 3
└── stimulus materials.pdf   # Visual stimuli used in the study
```

## Data Description

### 1. Fixation Data (`data_fixations.csv`)
Contains eye-tracking fixation data from the advertisement phase (Block 2). Used for Dirichlet regressions.
- Key variables: subject, trial, fix_type (Product/Face/Elsewhere), fix_time, Celebrity_Manipulation, Gaze_Manipulation

### 2. DDM Data (`data_ddm.csv`)
Contains data for logistic regressions and drift-diffusion model fitting.
- Key variables: subject, trial, value_left/right, response, reaction time, last_fixation, gaze_right, manipulations

### 3. Pupil Data (`data_pupil.csv`)
Contains processed pupil diameter measurements.
- Key variables: subject, trial, Celebrity_R, Time_event, pupil, Condition

## Stan Models

The repository includes several hierarchical drift-diffusion models implemented in Stan:
- Basic hierarchical DDMs (hDDM1-3)
- Attentional hierarchical DDMs (haDDM2-3)

These models were used to analyze the decision-making process and the influence of celebrity endorsements on consumer choices.

## Citation

If you use this data or code in your research, please cite:

```bibtex
@article{dambrogio2023celebrity,
  title={How celebrity status and gaze direction in ads drive visual attention to shape consumer decisions},
  author={D'Ambrogio, Pietro and Werksman, Nicolette and Platt, Michael L and Johnson, Eric J},
  journal={Psychology \& Marketing},
  volume={40},
  number={4},
  pages={741--755},
  year={2023},
  publisher={Wiley Online Library}
}
```

## License

Please refer to the original paper and authors for usage rights and permissions. 