# Theoretical Solar Cell Efficiencies Based on the Shockley Queisser limit
It calculates the theoretical solar cell parameters with options to change temperature, light intensity, and radiative efficiency, provides visualization tools. 

The calculation is based on the radiative limit or the Shockley Queisser limit

[Wikipedia: Shockley–Queisser limit] (https://en.wikipedia.org/wiki/Shockley–Queisser_limit)


#  Example Outputs
## Bandap vs Efficiencies of a single junction solar cell

Create a `SQlim` object, which calculates all the parameters and then call the method `plot('PCE')` to generate the bandgap vs efficiency plot. PCE: power conversion efficiency.
```python
SQ = SQlim()
SQ.plot('PCE')
```

<img src="/ExampleOutputFig/PCE.png" width="600">


The efficiencies will be in the `SQ.PCE` attribute, an numpy array. There's also a dictionary `SQ.paras` that stores all the characteristics as a key (string): characteristics (numpy array) pairs.


## Four important parameters in a single figure: 

Calling the method `plotall()` will generate a plot containing the 4 most important characteristics of a solar cell in the subplots
```python
SQ.plotall()
```
* Voc: Open-circuit voltage,
* Jsc: Short-circuit current density,
* FF: Fill factor

<img src="/ExampleOutputFig/ALL.png" width="800">

A method `get_paras(self, Eg, toPrint=True)` can be used to look up the results. For example, the following call would print the theoretical parameters for a 1.337 eV solar cell.
```python
SQ.get_paras(Eg = 1.337)
```
would print the following lines like these in the colsole:
```python
"""
Bandgap: 1.337 eV ; J0 = 2.64e-17 mA/cm^2

Voc = 1.079      V
Jsc = 35.14      mA/cm^2
FF  = 88.88      %
PCE = 33.703     %
"""
```



## Plot other characteristics

The `plot(para)` method can be used to generate different plots. Valid input `para` are `"Voc"`, `"Jsc"`, ``"FF"`, `"PCE"`, and `J0` (dark saturation current)
```python
SQ.plot('J0') # dark saturation current
SQ.plot('Voc') 
```

<img src="/ExampleOutputFig/J0.png" width="450"> <img src="/ExampleOutputFig/Voc.png" width="450">



## savedata
The data can be saved as a single .csv file
```python
SQ.saveall(savename = "SQ lim")
```

The data can be accessed here: [SQ limi data]("/SQ lim.csv")


#
# Visualize more interesting results

## Break down of the theoretical efficiency and the energy loss

The function provides in the script `E_loss`, which takes bandgap `Eg` as an input, and an optional input of an `SQlim` object, can be used to visualize the break down of energy loss. 

```python
E_loss(Eg = 1.337, SQ = SQ)
```

Shown here are the break down for a 1.337 eV solar cell, which has the maximum theoretical efficiency of 33.7 %.

<img src="/ExampleOutputFig/E_loss_1pt337eV.png" width="800">


##  Available Energies

The function E_available can be used to calculate and plot theoretical maximum available energies from a series of (mechanically stacked) solar cells with different Egs.

### Single-junction solar cell, 1.337 eV

```python
available_E(Egs = 1.337, SQ = SQ)
```

This is the similar to the one above but without the break down of energy loss.

<img src="/ExampleOutputFig/E_avail_1pt337eV.png" width="800">

  
    
      
## Multi-junction solar cells

The theoretical efficiencies of multijunction solar cells can be higher than single junction solar cells. The function `availableE` can actually take a list of different bandgaps and calculate the maximum possible efficiencies by using materials with these bandgaps. 

### Two bandgaps: 1.1 eV and 1.8 eV
```python
available_E(Egs = [1.1, 1.8], SQ = SQ)
```
The sum of the two sub-cells are higher than any single-junction solar cells.
<img src="/ExampleOutputFig/E_avail_2cells.png" width="800">

### Three bandgaps: 0.95 eV, 1.37 eV, 1.90 eV

```python
available_E(Egs = [0.95, 1.37, 1.90], SQ = SQ)
```
The sum of the efficiency are even higher.
<img src="/ExampleOutputFig/E_avail_3cells.png" width="800">


### Three bandgaps: Ge(0.65 eV), InGaAs (1.40 eV), InGaP (1.90 eV)

####  This bandgap-material combination is the example you can find on [Wikipedia's Multi-junction_solar_cell] (https://en.wikipedia.org/wiki/Multi-junction_solar_cell) page

<img src="/ExampleOutputFig/E_avail_3cells_InGaP_InGaAs_Ge.png" width="800">


#
#
#
# Different Conditions

The default conditions for calculating the theoretical limits are the standard conditions : Temperature `T = 300` K, 1-sun condition `intensity = 1.0`, and radiative efficiency `EQE_EL = 1.0` (100% external quantum efficiency for electroluminescence). 
```python
class SQlim(object):
    def __init__(self, T=300, EQE_EL=1.0, intensity=1.0):
        """
        T: temperature in K
        EQE_EL: radiative efficiency (EL quantum yield)
        intensity: light concentration, 1.0 = one Sun, 100 mW/cm^2
        """
```

We can calculate the efficiencies in different conditions by simply changing the input. 

Because of this flexibility, we can easily get an idea of how the change of these factors affect the theoretical efficiencies (or other characteristics fo interest).

## Different Temperature

#### The *theoretical* maximum possible efficiencies of solar cells could be higher at lower temperature.
The function `VaryTemp(T)` can do this calculation and plot the results.

```python
VaryTemp(T = [150, 200, 250, 300, 350, 400])
```

<img src="/ExampleOutputFig/VaryT.png" width="800">


## Different Light intensity (Solar concentrator)

####  The efficiencies are higher when the incident light intensity is higher. 
The function `VarySuns` does that caculation. 
```python
VarySuns(Suns = [1, 10, 100, 1000])
```

<img src="/ExampleOutputFig/VaryIntensity.png" width="800">


## Different radiative efficiency

###  The higher the EQE_EL (radiative efficiency), the higher the power conversion efficiency.

<img src="/ExampleOutputFig/VaryEQEEL.png" width="800">







