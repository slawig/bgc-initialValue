#!/usr/bin/env python
# -*- coding: utf8 -*

import itertools
import matplotlib.pyplot as plt
import numpy as np
import os

import metos3dutil.database.constants as DB_Constants
import metos3dutil.metos3d.constants as Metos3d_Constants
import metos3dutil.petsc.petscfile as petsc
from metos3dutil.plot.surfaceplot import SurfacePlot
import initialValue.constants as InitialValue_Constants
from initialValue.InitialValuePlot import InitialValuePlot


def main(orientation='lc2', fontsize=9, plotSpinup=False, plotNorm=False, plotRelationSpinupNorm=True, plotSurface=False):
    """
    Plot the results using different intial values for the spin up

    Parameters
    ----------
    orientation : str
        Orientation of the figure
    fontsize : int
        Fontsize used in the figure
    plotSpinup : bool, default: False
        If True, plot the figure of the spin up norm
    plotNorm : bool, default: False
        If True, plot the figures for the norm
    plotRelationSpinupNorm : bool, default: False
        If True, plot the figures of the relation between spin up norm and
        relative error
    """
    assert type(orientation) is str
    assert type(fontsize) is int and 0 < fontsize
    assert type(plotSpinup) is bool
    assert type(plotNorm) is bool
    assert type(plotRelationSpinupNorm) is bool

    initialValuePlot = InitialValuePlots(orientation=orientation, fontsize=fontsize)

    #Spin up norm plot
    if plotSpinup:
        parameterId = 0
        kwargs = {'N': {}, 'N-DOP': {}, 'NP-DOP': {}, 'NPZ-DOP': {}, 'NPZD-DOP': {}, 'MITgcm-PO4-DOP': {}}
        kwargs['N'] = {'orientation': 'lc2Legend', 'additionalSimulationIds': [(0, 'Reference'), (1806, 'Uniform'), (606, 'Normal'), (6, 'Lognormal'), (1206, 'OneBox')], 'legend_box': True, 'subplot_adjust': {'left': 0.195, 'bottom': 0.1775, 'right': 0.95, 'top': 0.8125}}
        kwargs['N-DOP'] = {'additionalSimulationIds': [(1, 'Reference'), (1906, 'Uniform'), (706, 'Normal'), (106, 'Lognormal'), (1306, 'OneBox'), (4806, 'Constant'), (4306, 'UniformRandom'), (3106, 'NormalRandom'), (2506, 'LognormalRandom'), (3706, 'OneBoxRandom')], 'legendLinestyles': True, 'legendLinestylesLoc': 'upper right'}
        kwargs['NP-DOP'] = {'orientation': 'lc2Legend', 'additionalSimulationIds': [(2, 'Reference'), (2006, 'Uniform'), (806, 'Normal'), (206, 'Lognormal'), (1406, 'OneBox'), (4905, 'Constant'), (4406, 'UniformRandom'), (3206, 'NormalRandom'), (2606, 'LognormalRandom'), (3806, 'OneBoxRandom')], 'legend_box': True, 'subplot_adjust': {'left': 0.195, 'bottom': 0.1775, 'right': 0.95, 'top': 0.8125}, 'legendLinestyles': True, 'legendLinestylesLoc': 'upper right'}
        kwargs['NPZ-DOP'] = {'additionalSimulationIds': [(3, 'Reference'), (2106, 'Uniform'), (906, 'Normal'), (306, 'Lognormal'), (1506, 'OneBox'), (5004, 'Constant'), (4506, 'UniformRandom'), (3306, 'NormalRandom'), (2706, 'LognormalRandom'), (3906, 'OneBoxRandom')]}
        kwargs['NPZD-DOP'] = {'additionalSimulationIds': [(4, 'Reference'), (2206, 'Uniform'), (1006, 'Normal'), (406, 'Lognormal'), (1606, 'OneBox'), (5103, 'Constant'), (4606, 'UniformRandom'), (3406, 'NormalRandom'), (2806, 'LognormalRandom'), (4006, 'OneBoxRandom')]}
        kwargs['MITgcm-PO4-DOP'] = {'additionalSimulationIds': [(5, 'Reference'), (2306, 'Uniform'), (1106, 'Normal'), (506, 'Lognormal'), (1706, 'OneBox'), (5202, 'Constant'), (4706, 'UniformRandom'), (3506, 'NormalRandom'), (2906, 'LognormalRandom'), (4106, 'OneBoxRandom')]}

        for metos3dModel in Metos3d_Constants.METOS3D_MODELS:
            initialValuePlot.plotSpinupData(metos3dModel, parameterId, **kwargs[metos3dModel])


    #Norm data plot
    if plotNorm:
        year = 10000
        normList = ['2']
        parameterIdList = [0]
        metos3dModelList = Metos3d_Constants.METOS3D_MODELS

        kwargs = {'N': {}, 'N-DOP': {}, 'NP-DOP': {}, 'NPZ-DOP': {}, 'NPZD-DOP': {}, 'MITgcm-PO4-DOP': {}}
        kwargs['N'] = {'orientation': 'lc2Legend', 'additionalSimulationIds': [(0, 'Reference'), (1806, 'Uniform'), (606, 'Normal'), (6, 'Lognormal'), (1206, 'OneBox')]}
        kwargs['N-DOP'] = {'additionalSimulationIds': [(1, 'Reference'), (1906, 'Uniform'), (706, 'Normal'), (106, 'Lognormal'), (1306, 'OneBox'), (4806, 'Constant'), (4306, 'UniformRandom'), (3106, 'NormalRandom'), (2506, 'LognormalRandom'), (3706, 'OneBoxRandom')]}
        kwargs['NP-DOP'] = {'orientation': 'lc2Legend', 'additionalSimulationIds': [(2, 'Reference'), (2006, 'Uniform'), (806, 'Normal'), (206, 'Lognormal'), (1406, 'OneBox'), (4905, 'Constant'), (4406, 'UniformRandom'), (3206, 'NormalRandom'), (2606, 'LognormalRandom'), (3806, 'OneBoxRandom')]}
        kwargs['NPZ-DOP'] = {'additionalSimulationIds': [(3, 'Reference'), (2106, 'Uniform'), (906, 'Normal'), (306, 'Lognormal'), (1506, 'OneBox'), (5004, 'Constant'), (4506, 'UniformRandom'), (3306, 'NormalRandom'), (2706, 'LognormalRandom'), (3906, 'OneBoxRandom')]}
        kwargs['NPZD-DOP'] = {'additionalSimulationIds': [(4, 'Reference'), (2206, 'Uniform'), (1006, 'Normal'), (406, 'Lognormal'), (1606, 'OneBox'), (5103, 'Constant'), (4606, 'UniformRandom'), (3406, 'NormalRandom'), (2806, 'LognormalRandom'), (4006, 'OneBoxRandom')]}
        kwargs['MITgcm-PO4-DOP'] = {'additionalSimulationIds': [(5, 'Reference'), (2306, 'Uniform'), (1106, 'Normal'), (506, 'Lognormal'), (1706, 'OneBox'), (5202, 'Constant'), (4706, 'UniformRandom'), (3506, 'NormalRandom'), (2906, 'LognormalRandom'), (4106, 'OneBoxRandom')]}
        for (norm, parameterId, metos3dModel) in list(itertools.product(normList, parameterIdList, metos3dModelList)):
            initialValuePlot.plotNormData(metos3dModel, parameterId, norm=norm, year=year, **kwargs[metos3dModel])


    #Plot relation between spin up norm and relative error for all parameter vectors
    if plotRelationSpinupNorm:
        year = 10000

        initialValues = [('vector', None, 'Uniform', 'set_mass', 'Uniform'), ('vector', None, 'Normal', 'set_mass', 'Normal'), ('vector', None, 'Lognormal', 'set_mass', 'Lognormal'), ('vector', None, 'OneBox', 'set_mass', 'OneBox'), ('constant', None, None, None, 'Constant'), ('vector', None, 'Uniform', 'random_mass', 'Uniform'), ('vector', None, 'Normal', 'random_mass', 'Normal'), ('vector', None, 'Lognormal', 'random_mass', 'Lognormal'), ('vector', None, 'OneBox', 'random_mass', 'OneBox')]
        zoom = {'N': [0, 1, 2], 'N-DOP': [0, 1, 2, 4, 5, 6, 7], 'NP-DOP': [0, 1, 2, 4, 5, 6, 7], 'NPZ-DOP': [0, 1, 2, 4, 5, 6, 7], 'NPZD-DOP': [0, 1, 2, 4, 5, 6, 7], 'MITgcm-PO4-DOP': [0, 1, 2, 4, 5, 6, 7]}

        for metos3dModel in Metos3d_Constants.METOS3D_MODELS:
            for filenamePrefix in ['', 'Zoom']:

                #Restrict the initial concentrations
                if filenamePrefix == 'Zoom':
                    initialValueParameter = [initialValues[i] for i in zoom[metos3dModel]]
                elif metos3dModel == 'N':
                    initialValueParameter = initialValues[:4]
                else:
                    initialValueParameter = initialValues

                kwargs = {'initialValueParameter': initialValueParameter, 'filenamePrefix': filenamePrefix}

                if filenamePrefix == 'Zoom':
                    kwargs.update({'xlabelpad': 12.0})

                #Log scale for x axis
                if filenamePrefix == '' and metos3dModel in ['NP-DOP', 'NPZ-DOP', 'NPZD-DOP']:
                    kwargs.update({'xscalelog': True})

                #Add a legend
                if metos3dModel in ['N', 'NP-DOP'] and filenamePrefix == '':
                    kwargs.update({'orientation': 'lc2Legend', 'legend_box': True})
                if metos3dModel in ['N-DOP', 'NP-DOP'] and filenamePrefix == '':
                    kwargs.update({'legendMarkers': True, 'legendMarkersLoc': 'lower left'})

                #Restriction of the data
                if filenamePrefix == 'Zoom' and metos3dModel == 'NP-DOP':
                    kwargs.update({'xrange': [1.392*10**(-5)+4.*10**(-9), 1.394*10**(-5)], 'yrange': [0.0, 10**(-3)]})
                if filenamePrefix == 'Zoom' and metos3dModel == 'NPZ-DOP':
                    kwargs.update({'xrange': [2.647*10**(-6), 2.649*10**(-6)], 'yrange': [0.0, 10**(-3)]})
                if filenamePrefix == 'Zoom' and metos3dModel == 'NPZD-DOP':
                    kwargs.update({'xrange': [8.4245*10**(-6), 8.4255*10**(-6)], 'yrange': [0.0, 10**(-3)]})

                #Subplot_adjust (if different from the standard values)
                if metos3dModel in ['N', 'NP-DOP'] and filenamePrefix == '':
                    kwargs.update({'subplot_adjust': {'left': 0.185, 'bottom': 0.185, 'right': 0.995, 'top': 0.8275}})
                elif metos3dModel in ['N-DOP'] and filenamePrefix == '':
                    kwargs.update({'subplot_adjust': {'left': 0.185, 'bottom': 0.225, 'right': 0.9375, 'top': 0.995}})
                elif metos3dModel in ['NPZD-DOP'] and filenamePrefix == '':
                    kwargs.update({'subplot_adjust': {'left': 0.185, 'bottom': 0.225, 'right': 0.995, 'top': 0.975}})
                elif filenamePrefix == 'Zoom' and metos3dModel in ['N']:
                    kwargs.update({'subplot_adjust': {'left': 0.185, 'bottom': 0.29, 'right': 0.995, 'top': 0.96}})
                elif filenamePrefix == 'Zoom':
                    kwargs.update({'subplot_adjust': {'left': 0.185, 'bottom': 0.29, 'right': 0.995, 'top': 0.995}})


                initialValuePlot.plotScatterSpinupNorm(metos3dModel, year=year, **kwargs)


    #Plot the tracer concentration difference at the surface
    if plotSurface:
        tracerCount = {'N': 'One_Tracer', 'N-DOP': 'Two_Tracer', 'NP-DOP': 'Three_Tracer', 'NPZ-DOP': 'Four_Tracer', 'NPZD-DOP': 'Five_Tracer', 'MITgcm-PO4-DOP': 'Two_Tracer'}
        surfacePlot = []
        surfacePlot += [('NP-DOP', 'Lognormal', 'random_mass', 0, False, False)]
        surfacePlot += [('NP-DOP', 'Normal', 'random_mass', 0, False, False)]
        surfacePlot += [('NPZ-DOP', 'OneBox', 'set_mass', 0, False, False)]
        #surfaceplot += [('NP-DOP', 'Lognormal', 'random_mass', 0, True, True), ('NP-DOP', 'Lognormal', 'random_mass', 0, True, False), ('NP-DOP', 'Lognormal', 'random_mass', 0, False, False)]
        #surfacePlot += [('NPZ-DOP', 'Lognormal', 'random_mass', 0, True, True)] #, ('NPZ-DOP', 'Lognormal', 'random_mass', 0, True, False), ('NPZ-DOP', 'Lognormal', 'random_mass', 0, False, False)]
        #surfacePlot += [('NPZ-DOP', 'Normal', 'random_mass', 0, True, True)] #, ('NPZ-DOP', 'Normal', 'random_mass', 0, True, False), ('NPZ-DOP', 'Normal', 'random_mass', 0, False, False)]
        #surfacePlot += [('NPZ-DOP', 'OneBox', 'random_mass', 0, True, True)] #, ('NPZ-DOP', 'OneBox', 'random_mass', 0, True, False), ('NPZ-DOP', 'OneBox', 'random_mass', 0, False, False)]
        #surfacePlot += [('NPZ-DOP', 'Uniform', 'random_mass', 0, True, True)] #, ('NPZ-DOP', 'Uniform', 'random_mass', 0, True, False), ('NPZ-DOP', 'Uniform', 'random_mass', 0, False, False)]
        #surfacePlot += [('NPZ-DOP', 'OneBox', 'set_mass', 0, True, True)] #, ('NPZ-DOP', 'OneBox', 'set_mass', 0, True, False), ('NPZ-DOP', 'OneBox', 'set_mass', 0, False, False)]
        #surfacePlot += [('NPZD-DOP', 'Lognormal', 'random_mass', 0, True, True)] #, ('NPZ-DOP', 'Lognormal', 'random_mass', 0, True, False), ('NPZ-DOP', 'Lognormal', 'random_mass', 0, False, False)]
        #surfacePlot += [('NPZD-DOP', 'Normal', 'random_mass', 0, True, True)] #, ('NPZ-DOP', 'Normal', 'random_mass', 0, True, False), ('NPZ-DOP', 'Normal', 'random_mass', 0, False, False)]
        #surfacePlot += [('NPZD-DOP', 'OneBox', 'random_mass', 0, True, True)] #, ('NPZ-DOP', 'OneBox', 'random_mass', 0, True, False), ('NPZ-DOP', 'OneBox', 'random_mass', 0, False, False)]
        #surfacePlot += [('NPZD-DOP', 'Uniform', 'random_mass', 0, True, True)] #, ('NPZ-DOP', 'Uniform', 'random_mass', 0, True, False), ('NPZ-DOP', 'Uniform', 'random_mass', 0, False, False)]
        #surfacePlot += [('NPZD-DOP', 'OneBox', 'set_mass', 0, True, True)] #, ('NPZ-DOP', 'OneBox', 'set_mass', 0, True, False), ('NPZ-DOP', 'OneBox', 'set_mass', 0, False, False)]

        for (model, distribution, tracerDistribution, tracerNum, tracerDifference, relativeError) in surfacePlot:
            filenameFigure = os.path.join(InitialValue_Constants.PATH_FIGURE, 'Spinup', model, InitialValue_Constants.PATTERN_FIGURE_SURFACE.format(model, 1, 0, distribution, tracerDistribution, '{}', '{}', relativeError, tracerDifference))
            filenameTracer = os.path.join(InitialValue_Constants.PATH, 'Data', model, InitialValue_Constants.PATH_DISTRIBUTION.format(distribution), InitialValue_Constants.PATH_TRACER_DISTRIBUTION.format(tracerDistribution), tracerCount[model], 'InitialTracer_{:0>3d}'.format(tracerNum), '1dt', 'Tracer', Metos3d_Constants.PATTERN_TRACER_OUTPUT)
            filenameTracerReference = os.path.join(InitialValue_Constants.PATH, 'Data', model, 'Parameter_{:0>3d}'.format(0), InitialValue_Constants.PATH_CONCENTRATION_TYP['constant'], InitialValue_Constants.PATH_DIFFERENT_TRACER[len(Metos3d_Constants.METOS3D_MODEL_TRACER[model])], 'InitialTracer_{:0>3d}'.format(0), 'Tracer', Metos3d_Constants.PATTERN_TRACER_OUTPUT)

            initialValuePlot.plotTracerConcentrationSurface(model, filenameFigure, filenameTracer, filenameTracerReference=filenameTracerReference, tracerDifference=tracerDifference, relativeError=relativeError, plotSlice=False, slicenum=None)



class InitialValuePlots():
    """
    Preparation of the plots for the results using different initial values

    Attributes
    ----------
    orientation : str
        Orientation of the figure
    fontsize : int
        Fontsize used in the figure
    """

    def __init__(self, orientation='lc2', fontsize=9):
        """
        Constructs the environment to plot the data using different initial
        values for the spin up calculation.

        Parameter
    ----------
        orientation : str, default: lc2
            Orientation of the figure
        fontsize : int, default: 9
            Fontsize used in the figure
        """
        assert type(orientation) is str
        assert type(fontsize) is int and 0 < fontsize

        self._orientation = orientation
        self._fontsize = fontsize

        self.__initialValuePlot = InitialValuePlot(orientation=self._orientation, fontsize=self._fontsize)


    def closeDatabaseConnection(self):
        """
        Close the connection of the database
        """
        self.__initialValuePlot.closeDatabaseConnection()


    def plotSpinupData(self, metos3dModel, parameterId, ncol=3, subPlot=False, **kwargs):
        """
        Plot the spin up norm.

        Plot the spin up norm for the given biogeochemical model and
        parameterId. The plot includes the spin up norm for the reference
        simulation (using the standard constant initial value) and for the
        additional simulations defined in additionalSimulationIds of kwargs.

        Parameters
        ----------
        metos3dModel : str
            Name of the biogeochemical model
        parameterId : int
            Id of the parameter of the latin hypercube example
        ncol : int, default: 3
            Number of columns for the legend
        subPlot : bool, default: True
            If the value is True, an enlargment of the last 2000 model years
            of the spin up norm is inserted as a extra subplot. For the value
            False, no subplot is added.
        **kwargs : dict
            Additional keyword arguments with keys:

            orientation : str, optional
                Orientation of the figure
            axesResultSmall : list [float], optional
                Dimensions of the subplot
            subPlotModelYear : int, optional
                Starting model year of the sub plot
            subplot_adjust : dict [str, float], optional
                Adjustment of the subplot using the keys left, bottom, right
                and top
            legend_box : bool, optional
                If the value is True, plot the legend of the plot using an bbox
                above the plot
            additionalSimulationIds : list [tuple], optional
                List for additional spin up plots using the simulationId and
                label defined in the tuples
            filenamePrefix : str or None, optional
                Prefix of the filename for the figure

        NOTES
        -----
        The figure is saved in the directory defined in
        initialValue.constants.PATH_FIGURE.
        """
        assert metos3dModel in Metos3d_Constants.METOS3D_MODELS
        assert parameterId in range(0, InitialValue_Constants.PARAMETERID_MAX+1)
        assert type(ncol) is int and 0 < ncol
        assert type(subPlot) is bool

        #Parse keyword arguments
        orientation = kwargs['orientation'] if 'orientation' in kwargs and type(kwargs['orientation']) is str else self._orientation
        axesResultSmall = kwargs['axesResultSmall'] if 'axesResultSmall' in kwargs and type(kwargs['axesResultSmall']) is list and len(kwargs['axesResultSmall']) == 4 else [.69, .5, .2, .2]
        subPlotModelYear = kwargs['subPlotModelYear'] if 'subPlotModelYear' in kwargs and type(kwargs['subPlotModelYear']) is int and 0 < kwargs['subPlotModelYear'] else 8000
        subplot_adjust = kwargs['subplot_adjust'] if 'subplot_adjust' in kwargs and type(kwargs['subplot_adjust']) is dict and 'left' in kwargs['subplot_adjust'] and 'bottom' in kwargs['subplot_adjust'] and 'right' in kwargs['subplot_adjust'] and 'top' in kwargs['subplot_adjust'] else {'left': 0.195, 'bottom': 0.2175, 'right': 0.95, 'top': 0.995}
        additionalSimulationIds = kwargs['additionalSimulationIds'] if 'additionalSimulationIds' in kwargs and type(kwargs['additionalSimulationIds'] is list) else []
        filenameSpinup = os.path.join(InitialValue_Constants.PATH_FIGURE, 'Spinup', metos3dModel, kwargs['filenamePrefix'] + InitialValue_Constants.PATTERN_FIGURE_SPINUP.format(metos3dModel, parameterId) if 'filenamePrefix' in kwargs and type(kwargs['filenamePrefix']) is str else InitialValue_Constants.PATTERN_FIGURE_SPINUP.format(metos3dModel, parameterId))
        kwargsLegend = {'legend': kwargs['legend']} if 'legend' in kwargs and type(kwargs['legend']) is bool else {}
        kwargsLegend.update({'legendLinestyles': kwargs['legendLinestyles']} if 'legendLinestyles' in kwargs and type(kwargs['legendLinestyles']) is bool else {})
        kwargsLegend.update({'legendLinestylesLoc': kwargs['legendLinestylesLoc']} if 'legendLinestylesLoc' in kwargs and type(kwargs['legendLinestylesLoc']) is str else {})
        kwargsLegend.update({'labelspacing': 0.3, 'handlelength': 0.8, 'handletextpad': 0.4, 'columnspacing': 1.0, 'borderaxespad': 0.3, 'borderpad': 0.2})

        #Create spin up norm plot
        self.__initialValuePlot._init_plot(orientation=orientation, fontsize=self._fontsize)
        self.__initialValuePlot.plot_spinup_data(ncol=ncol, simulationIds=additionalSimulationIds, subPlot=subPlot, axesResultSmall=axesResultSmall, subPlotModelYear=subPlotModelYear, **kwargsLegend)

        if 'legend_box' in kwargs and type(kwargs['legend_box']) is bool and kwargs['legend_box']:
            self.__initialValuePlot.set_legend_box(bbox_to_anchor=(-0.08, 1.02, 1.15, 0.2), ncol=ncol, labelspacing=kwargsLegend['labelspacing'], handlelength=kwargsLegend['handlelength'], handletextpad=kwargsLegend['handletextpad'], columnspacing=kwargsLegend['columnspacing'], borderaxespad=kwargsLegend['borderaxespad'], borderpad=kwargsLegend['borderpad'])

        self.__initialValuePlot.set_subplot_adjust(left=subplot_adjust['left'], bottom=subplot_adjust['bottom'], right=subplot_adjust['right'], top=subplot_adjust['top'])

        self.__initialValuePlot.savefig(filenameSpinup)
        self.__initialValuePlot.close_fig()



    def plotNormData(self, model, parameterId, norm='2', trajectory='', year=None, ncol=3, **kwargs):
        """
        Plot the relative error

        Plot the relative error of the spin up for the given biogeochemical
        model and parameterId. The plot includes the relative error for the
        simulations defined in additionalSimulationIds of kwargs.

        Parameters
        ----------
        model : str
            Name of the biogeochemical model
        parameterId : int
            Id of the parameter of the latin hypercube example
        norm : str, default: 2
            Descriptive string for the norm to be used
            (see util.metos3dutil.database.constants.NORM)
        trajectory : str, default: ''
            Use for '' the norm only at the first time point in a model year
            and use for 'trajectory' the norm over the whole trajectory
        year : int or None, default: None
            Model year used for the reference solution
            If None, use the same year for the reference solution,
            respectively
        ncol : int, default: 3
            Number of columns for the legend
        **kwargs : dict
            Additional keyword arguments with keys:

            subplot_adjust : dict [str, float], optional
                Adjustment of the subplot using the keys left, bottom, right
                and top
            legend_box : bool, optional
                If the value is True, plot the legend of the plot using an bbox
                above the plot
            additionalSimulationIds : list [tuple [int]], optional
                List for additional spin up plots using the simulationId and
                timestep defined in the tuples
            filenamePrefix : str or None, optional
                Prefix of the filename for the figure

        NOTES
        -----
        The figure is saved in the directory defined in
        initialValue.constants.PATH_FIGURE.
        """
        assert model in Metos3d_Constants.METOS3D_MODELS
        assert parameterId in range(0, InitialValue_Constants.PARAMETERID_MAX+1)
        assert trajectory in ['', 'trajectory']
        assert norm in DB_Constants.NORM
        assert year is None or type(year) is int and 0 <= year
        assert type(ncol) is int and 0 < ncol

        #Parse keyword arguments
        subplot_adjust = kwargs['subplot_adjust'] if 'subplot_adjust' in kwargs and type(kwargs['subplot_adjust']) is dict and 'left' in kwargs['subplot_adjust'] and 'bottom' in kwargs['subplot_adjust'] and 'right' in kwargs['subplot_adjust'] and 'top' in kwargs['subplot_adjust'] else {'left': 0.185, 'bottom': 0.2175, 'right': 0.95, 'top': 0.995}
        additionalSimulationIds = kwargs['additionalSimulationIds'] if 'additionalSimulationIds' in kwargs and type(kwargs['additionalSimulationIds'] is list) else []
        filenameNorm = os.path.join(InitialValue_Constants.PATH_FIGURE, 'Norm', kwargs['filenamePrefix'] + InitialValue_Constants.PATTERN_FIGURE_NORM.format(trajectory, norm, model, parameterId) if 'filenamePrefix' in kwargs and type(kwargs['filenamePrefix']) is str else InitialValue_Constants.PATTERN_FIGURE_NORM.format(trajectory, norm, model, parameterId))
        kwargsLegend = {'handlelength': 1.0, 'handletextpad': 0.4, 'columnspacing': 1.0, 'borderaxespad': 0.3}

        #Create scatter plot
        self.__initialValuePlot._init_plot(orientation=self._orientation, fontsize=self._fontsize)
        self.__initialValuePlot.plot_tracer_norm_data(additionalSimulationIds[0][0], additionalSimulationIds[1:], ncol=ncol, year=year, **kwargsLegend)

        if 'legend_box' in kwargs and type(kwargs['legend_box']) is bool and kwargs['legend_box']:
            self.__initialValuePlot.set_legend_box(ncol=ncol)

        self.__initialValuePlot.set_subplot_adjust(left=subplot_adjust['left'], bottom=subplot_adjust['bottom'], right=subplot_adjust['right'], top=subplot_adjust['top'])
        self.__initialValuePlot.savefig(filenameNorm)
        self.__initialValuePlot.close_fig()


    def plotScatterSpinupNorm(self, model, year=10000, norm='2', trajectory='', **kwargs):
        """
        Plot the spin up norm against the norm

        Plot the spin up norm value against the norm value for the given model
        year using the different initial values

        Parameters
        ----------
        model : str
            Name of the biogeochemical model
        year : int, default: 10000
            Used model year of the spin up (for the spin up is the previous
            model year used)
        norm : str, default: 2
            Descriptive string for the norm to be used
            (see util.metos3dutil.database.constants.NORM)
        trajectory : str, default: ''
            Use for '' the norm only at the first time point in a model year
            and use for 'trajectory' the norm over the whole trajectory
        **kwargs : dict
            Additional keyword arguments with keys:

            orientation : str, optional
                Orientation of the figure
            initialValueParameter: list [tuple]
                List with tuples of configuration for different initial values
                consisting of (concentrationTyp, differentTracer, distribution,
                tracerDistribution, label)
            legend_box : bool
                If the value is True, plot the legend of the plot using an bbox
                above the plot
            ncol : int, default: 3
                Number of columns for the legend
            filenamePrefix : str or None, optional
                Prefix of the filename for the figure
            legendMarkers : bool
                If True, create legend with the different markers defining the
                different initial concentrations
            legendMarkersLoc : str
                Location of the legend for the different markers
            xscalelog : bool
                If True, use a logarithmic scale for the x axis
            xlabelpad : float
                Spacing in points from the axes bounding box including ticks
                and tick labels.
            xrange : list [float]
                Restriction of the spin-up value to the given range. If a value
                is None, is range is open on this side
            yrange : list [float]
                Restriction of the relative error value to the given range. If
                a value is None, is range is open on this side

        NOTES
        -----
        The figure is saved in the directory defined in
        initialValue.constants.PATH_FIGURE.
        """
        assert model in Metos3d_Constants.METOS3D_MODELS
        assert type(year) is int and 0 <= year
        assert norm in DB_Constants.NORM
        assert trajectory in ['', 'Trajectory']

        #Parse keyword arguments
        orientation = kwargs['orientation'] if 'orientation' in kwargs and type(kwargs['orientation']) is str else self._orientation
        initialValueParameter = kwargs['initialValueParameter'] if 'initialValueParameter' in kwargs and type(kwargs['initialValueParameter']) is list else [('constant', len(Metos3d_Constants.METOS3D_MODEL_TRACER[model]), None, None, 'Constant'), ('vector', len(Metos3d_Constants.METOS3D_MODEL_TRACER[model]), 'Lognormal', 'set_mass', 'Lognormal'), ('vector', len(Metos3d_Constants.METOS3D_MODEL_TRACER[model]), 'Normal', 'set_mass', 'Normal'), ('vector', len(Metos3d_Constants.METOS3D_MODEL_TRACER[model]), 'OneBox', 'set_mass', 'OneBox'), ('vector', len(Metos3d_Constants.METOS3D_MODEL_TRACER[model]), 'Uniform', 'set_mass', 'Uniform')]
        subplot_adjust = kwargs['subplot_adjust'] if 'subplot_adjust' in kwargs and type(kwargs['subplot_adjust']) is dict and 'left' in kwargs['subplot_adjust'] and 'bottom' in kwargs['subplot_adjust'] and 'right' in kwargs['subplot_adjust'] and 'top' in kwargs['subplot_adjust'] else {'left': 0.185, 'bottom': 0.225, 'right': 0.995, 'top': 0.995}
        ncol = kwargs['ncol'] if 'ncol' in kwargs and type(kwargs['ncol']) is int and kwargs['ncol'] > 0 else 3
        ncol = kwargs['ncol'] if 'ncol' in kwargs and type(kwargs['ncol']) is int and kwargs['ncol'] > 0 else 3
        handlelength = kwargs['handlelength'] if 'handlelength' in kwargs and type(kwargs['handlelength']) is float and 0.0 < kwargs['handlelength'] else 0.5
        handletextpad = kwargs['handlepadtext'] if 'handlepadtext' in kwargs and type(kwargs['handlepadtext']) is float and 0.0 < kwargs['handlepadtext'] else 0.4
        kwargsLegend = {'legendMarkers': kwargs['legendMarkers']} if 'legendMarkers' in kwargs and type(kwargs['legendMarkers']) is bool else {}
        kwargsLegend.update({'legendMarkersLoc': kwargs['legendMarkersLoc'] if 'legendMarkersLoc' in kwargs and type(kwargs['legendMarkersLoc']) is str else 'best'})
        kwargsLegend.update({'xscalelog': kwargs['xscalelog'] if 'xscalelog' in kwargs and type(kwargs['xscalelog']) is bool else False})
        kwargsLegend.update({'xlabelpad': kwargs['xlabelpad'] if 'xlabelpad' in kwargs and type(kwargs['xlabelpad']) is float else None})
        kwargsLegend.update({'xrange': kwargs['xrange']} if 'xrange' in kwargs and type(kwargs['xrange']) is list else {})
        kwargsLegend.update({'yrange': kwargs['yrange']} if 'yrange' in kwargs and type(kwargs['yrange']) is list else {})

        filename = os.path.join(InitialValue_Constants.PATH_FIGURE, 'Norm', kwargs['filenamePrefix'] + InitialValue_Constants.PATTERN_FIGURE_SPINUP_NORM.format(trajectory, norm, model) if 'filenamePrefix' in kwargs and type(kwargs['filenamePrefix']) is str else InitialValue_Constants.PATTERN_FIGURE_SPINUP_NORM.format(trajectory, norm, model))

        #Create scatter plot
        self.__initialValuePlot._init_plot(orientation=orientation, fontsize=self._fontsize)
        self.__initialValuePlot.plot_scatter_spinup_norm(model, initialValueParameter=initialValueParameter, year=year, norm=norm, trajectory=trajectory, **kwargsLegend)

        if 'legend_box' in kwargs and type(kwargs['legend_box']) is bool and kwargs['legend_box']:
            self.__initialValuePlot.set_legend_box(ncol=ncol, handlelength=handlelength, handletextpad=handletextpad)

        self.__initialValuePlot.set_subplot_adjust(left=subplot_adjust['left'], bottom=subplot_adjust['bottom'], right=subplot_adjust['right'], top=subplot_adjust['top'])

        self.__initialValuePlot.savefig(filename)
        self.__initialValuePlot.close_fig()


    def plotTracerConcentrationSurface(self, model, filenameSurface, filenameTracer, filenameTracerReference=None, tracerDifference=False, relativeError=False, cmap=None, plotSurface=True, plotSlice=False, slicenum=None, orientation='etnasp'):
        """
        Plot the tracer concentration for the given layer

        Parameters
        ----------
        model : str
            Name of the biogeochemical model
        filenameSurface : str
            Filename of the figure
        filenameTracer : str
            Pattern of the tracer filename including the path
        filenameTracerReference : str or None, default: None
            Pattern of the tracer filename used as reference tracer
            including the path
        tracerDifference : bool, default: False
            If True, plot the tracer concentration difference between the
            concentration calculated with the given time step and the
            concentration calculated with the time step 1dt
        relativeError : bool, default: False
            If True
        cmap :

        plotSurface : bool, default: True
            If True, plot the tracer concentration at the surface
        plotSlice : bool, default: True
            If True, plot the slices of the tracer concentration
        slicenum : None or list [int], default: None
            The list slice of the tracer concentration for the given boxes.
        orientation : str, default: 'etnasp'
            Orientation of the figure

        NOTES
        -----
        The figure is saved in the directory defined in
        timesteps.constants.PATH_FIGURE.
        """
        assert type(filenameTracer) is str
        assert type(filenameSurface) is str
        assert model in Metos3d_Constants.METOS3D_MODELS
        assert type(tracerDifference) is bool
        assert type(relativeError) is bool
        assert type(plotSurface) is bool
        assert type(plotSlice) is bool
        assert slicenum is None or isinstance(slicenum, list)
        assert filenameTracerReference is None and not tracerDifference and not relativeError or type(filenameTracerReference) is str

        #Check if the tracer exists
        tracerExists = True
        for tracer in Metos3d_Constants.METOS3D_MODEL_TRACER[model]:
            tracerExists = tracerExists and os.path.exists(filenameTracer.format(tracer)) and os.path.isfile(filenameTracer.format(tracer))

        if tracerExists:
            #Read tracer concentration
            tracerConcentration = np.zeros(shape=(Metos3d_Constants.METOS3D_VECTOR_LEN, len(Metos3d_Constants.METOS3D_MODEL_TRACER[model])))
            i = 0
            for tracer in Metos3d_Constants.METOS3D_MODEL_TRACER[model]:
                tracerConcentration[:,i] = petsc.readPetscFile(filenameTracer.format(tracer))
                i += 1
        else:
            #Missing tracer file
            assert False

        #Calculate the norm of the tracer concentration vector
        if tracerDifference or relativeError:
            tracer1dt = np.zeros(shape=(Metos3d_Constants.METOS3D_VECTOR_LEN, len(Metos3d_Constants.METOS3D_MODEL_TRACER[model])))
            i = 0
            for tracer in Metos3d_Constants.METOS3D_MODEL_TRACER[model]:
                tracer1dt[:,i] = petsc.readPetscFile(filenameTracerReference.format(tracer))
                i += 1

        if relativeError:
            normValue = np.linalg.norm(tracer1dt)
            print('Model: {}\nNorm rel: {}'.format(model, np.linalg.norm(tracerConcentration - tracer1dt) / normValue))
        else:
            normValue = 1.0

        #Plot the tracer concentration for every tracer
        i = 0
        for tracer in Metos3d_Constants.METOS3D_MODEL_TRACER[model]:
            if tracerDifference:
                v1d = np.divide(np.fabs(tracerConcentration[:,i] - tracer1dt[:,i]), normValue)
            else:
                v1d = np.divide(tracerConcentration[:,i], normValue)

            print('Tracer: {} Difference: {} rel: {} Norm: {}'.format(tracer, tracerDifference, relativeError, np.linalg.norm(v1d)))

            for depth in range(15):
                surface = SurfacePlot(orientation=orientation)
                surface.init_subplot(1, 2, orientation=orientation, gridspec_kw={'width_ratios': [9,5]})

                #Plot the surface concentration
                meridians = None if slicenum is None else [np.mod(Metos3d_Constants.METOS3D_GRID_LONGITUDE * x, 360) for x in slicenum]
                #cntr = surface.plot_surface(v1d, projection='robin', levels=50, ticks=plt.LogLocator(), format='%.1e', extend='both', meridians=meridians, colorbar=True)
                cntr = surface.plot_surface(v1d, projection='robin', levels=50, depth=depth, ticks=plt.LinearLocator(6), format='%.1e', extend='max', meridians=meridians, colorbar=True)
                #cntr = surface.plot_surface(v1d, projection='robin', levels=50, ticks=plt.LinearLocator(6), format='%.1e', pad=0.05, extend='max', meridians=meridians, colorbar=False)

                #Plot the slice plan of the concentration
                if plotSlice and slicenum is not None:
                    surface.set_subplot(0,1)
                    for s in slicenum:
                        surface.plot_slice(v1d, s, levels=50, ticks=plt.LinearLocator(6), format='%.1e', pad=0.02, extend='max', colorbar=False)

                plt.tight_layout(pad=0.05, w_pad=0.15)
                #cbar = surface._fig.colorbar(cntr, ax=surface._axes[0], format='%.1e', ticks=plt.LinearLocator(5), pad=0.02, aspect=40, extend='max', orientation='horizontal', shrink=0.8)

                surface.savefig(filenameSurface.format(tracer, depth))
                plt.close('all')
            i = i + 1



if __name__ == '__main__':
    main()

