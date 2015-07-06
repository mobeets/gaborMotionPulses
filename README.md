## Quickstart

1) Clone [gaborMotionPulses](https://github.com/mobeets/gaborMotionPulses) and [mASD](https://github.com/mobeets/mASD) to your home folder.

```
$ cd ~
$ git clone https://github.com/mobeets/gaborMotionPulses.git
$ git clone https://github.com/mobeets/mASD.git
```

2) Download [cbrewer](http://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab) to your home folder.

3) Add `nancyNeuronFiles, nancyStimFiles, patNeuronFiles, patStimFiles` to your home folder.

4) Open up Matlab and add paths.

```
>> cd ~/gaborMotionPulses
>> addpath ../mASD
>> addpath ../cbrewer
```

5) Fit ASD, ML, and Flat weights for all cells and behavior.

```
>> fitdir = 'runname';
>> fitAllSTRFs(fitdir, false, 'behavior ASD Flat ML');
>> fitAllSTRFs(fitdir, true, 'behavior ASD Flat ML');
>> fitAllSTRFs(fitdir, false, 'cells ASD Flat ML');
>> fitAllSTRFs(fitdir, true, 'cells ASD Flat ML');
```

The results now live in `~/gaborMotionPulses/data/runname-pat` and `~/gaborMotionPulses/data/runname-nancy`.

6) Load all ASD fits into struct array.

```
>> fitdir = 'runname';
>> vn = tools.makeFitSummaries(['data/' fitdir '-nancy/fits'], true, 'ASD');
>> vp = tools.makeFitSummaries(['data/' fitdir '-pat/fits'], false, 'ASD');
>> vs = [vp vn];
```
