## Macros to create plots used in the poster for PET is Wonderful 2023:

To get the count rate vs. activity plots run:
[TracerStudy.py](./TracerStudy.py)

To get the csv files and plot count rate vs. activity for different crystal materials on one figure run:
[CompareRateVsActivity.py](./CompareRateVsActivity.py)
**IMPORTANT:** The code assumes that the hits.csv files corresponding to different materials are already created using [TracerStudy.py](./TracerStudy.py) and placed in different directories so it may not work out of the box and instead you will have to update the paths

To get the poster-like nice plots using csv files created with the file above run: 
[PlotRateComparison.py](./PlotRateComparison.py)

To get the NECR vs. Length plots run:
[NECRvsLength.py](./NECRvsLength.py) 

To get the plot showing NECR for different tracers but with the same crystal material run: 
[PlotTracerNECR.py](./PlotTracerNECR.py)

To get nice poster-style plots comparing NECR vs. Length for different crystal materials run:
[CompareNECRvsLength.py](./CompareNECRvsLength.py)
**IMPORTANT:** The code assumes that the hits.csv files corresponding to different materials are already created using [NECRvsLength.py](./NECRvsLength.py) and placed in different directories so it may not work out of the box and instead you will have to update the paths




