Read Me for:
D'Ambrogio, Werksman, Platt, and Johnson 2021

———————————————————————————————Data——————————————————————————————————————————
File: data_fixation (included as .csv)  
Contains the fixation data of block 2 (advertisement phase). We used this dataset to fit the dirichlet regressions.
Each row is a fixation type.

subject	    	    = subject ID
trial			= trial number for subject
fix_type               = type of fixation (Product, Face, Elsewhere)
fix_time               = proportional fixation time
Celebrity_Manipulation = type of celebrity advertisement (Celebrity, Non-Celebrity, Snack Only) of the option on the right
Gaze_Manipulation      = type of gaze-cueing advertisement (Gaze on Product, Gaze on Viewer, Snack Only) of the option on the right
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
File: data_ddm (included as .csv)  
Contains the data used to run the logistic regressions and to fit the drift-diffusion model.
Each row is a trial.

subject		= subject ID
sbj_gender             = gender
trial			= trial number for subject
value_left		= subject’s subjective value for snack on the left-hand side of the screen
value_right		= subject’s subjective value for snack on the right-hand side of the screen
RV			= relative value, computed as the subjective value difference between the snack on the right and snack on the left.
response		= subject’s response for the current trial (left = 0, right = not 1)
rt			= reaction time in the current trial (milliseconds)
last_fixation		= subject’s direction of his last fixation in the current trial (Right item, Left item)
gaze_right		= proportional fixation duration on the right snack, computed by summing the fixation duration on the right item and dividing by the total trial duration
Celebrity_Manipulation = type of celebrity advertisement (Celebrity, Non-Celebrity, Snack Only) of the option on the right
Gaze_Manipulation      = type of gaze-cueing advertisement (Gaze on Product, Gaze on Viewer, Snack Only) of the option on the right
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
File: data_pupil (included as .csv)
Contains processed pupil data. Each row is a time step.

subject		= subject ID
trial			= trial number
Celebrity_R            = type of celebrity advertisement (Celebrity, Non-Celebrity, Snack Only) of the option on the right
Time_event             = time step in milliseconds
pupil                  = baseline-corrected pupil diameter, computed subtracting the mean value of the 200ms before and after the decision from each value within the analysis time window.
Condition              = type of celebrity advertisement (Celebrity, Non-C) of the option on the right. The level "Non-C" is the union of Non-Celebrity and Snack Only.
- - - - — - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

