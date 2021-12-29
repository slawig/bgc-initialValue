#!/usr/bin/env python
# -*- coding: utf8 -*

import matplotlib.pyplot as plt
import os
import numpy as np

import metos3dutil.database.constants as DB_Constants
import metos3dutil.metos3d.constants as Metos3d_Constants
from metos3dutil.plot.plot import Plot
import initialValue.constants as InitialValue_Constants
from initialValue.InitialValueDatabase import InitialValue_Database


class InitialValuePlot(Plot):
    """
    Creation of plots using different initial values
    """

    def __init__(self, orientation='lc2', fontsize=9, dbpath=InitialValue_Constants.DB_PATH, cmap=None, completeTable=True):
        """
        Constructs the environment to plot the data using different initial values

        Parameter
        ----------
        orientation : str, default: lc2
            Orientation of the figure
        fontsize : int, default: 9
            Fontsize used in the figure
        dbpath : str, default: initialValue.constants.DB_PATH
            Path to the sqlite database
        cmap : matplotlib.colors.Colormap or None, default: None
            Colormap used in the surface plot to visualize tracer
            concentrations
        completeTable : bool, default: True
            If the value is True, use all columns (even columns with value
            None) in SELECT queries on the database
        """
        assert type(orientation) is str
        assert type(fontsize) is int and 0 < fontsize
        assert os.path.exists(dbpath) and os.path.isfile(dbpath)
        assert type(completeTable) is bool

        Plot.__init__(self, cmap=cmap, orientation=orientation, fontsize=fontsize)
        self._database = InitialValue_Database(dbpath=dbpath, completeTable=completeTable)

        self._linestyle = ['solid', (0, (5, 10))]
        self._linestyleLabel = ['Default', 'Random']
        self._marker = ['.', '1']


    def closeDatabaseConnection(self):
        """
        Close the connection of the database
        """
        self._database.close_connection()


    def _linestyle_color(self, label):
        """
        Mapping the label to a linestyle and color

        Parameters
        ----------
        label : str
            Label of the data

        Returns
        -------
        tuple [int, bool]
            Tuple of the linestyle index, color index and a flag to plot the
            legend of the linestyle
        """
        assert type(label) is str

        legend = True

        #Set linestyle
        linestyle = 0
        if label.endswith('Random') or label == 'Constant':
            linestyle = 1
            legend = False

        #Set color
        colorIndex = -1
        if label == 'Constant':
            colorIndex = 4
        elif label.startswith('Uniform'):
            colorIndex = 1
        elif label.startswith('Normal'):
            colorIndex = 0
        elif label.startswith('Lognormal'):
            colorIndex = 2
        elif label.startswith('OneBox'):
            colorIndex = 3

        return (linestyle, colorIndex, legend)


    def plot_spinup_data(self, ncol=3, simulationIds=[], subPlot=False, axesResultSmall=[.61, .30, .3, .34], subPlotModelYear=8000, **kwargs):
        """
        Plot the spinup for the given simulationIds

        Parameter
        ---------
        ncol : int, default: 3
            Number of columns for the legend
        simulationIds : list [tuple], default: []
            List for additional spin up plots using the simulationId and label
            defined in the tuples. Add the label to to legend if the label is
            not ''
        subPlot : bool, default: False
            If the value is True, an enlargment of the last 2000 model years
            of the spin up norm is inserted as a extra subplot. For the value
            False, no subplot is added
        axesResultSmall : list [float], default: [.61, .30, .3, .34]
            Dimensions of the subplot
        subPlotModelYear : int, default: 8000
            Start model year for the subplot
        **kwargs : dict
            Additional keyword arguments with keys:

            legend : bool
                If True, plot the legend
            handlelength : float
                The length of the legend handles, in font-size units.
            handletextpad : float
                The pad between the legend handle and text, in font-size units.
            columnspacing : float
                The spacing between columns, in font-size units.
            borderaxespad : float
                The pad between the axes and legend border, in font-size units.
            legendLinestyles : bool
                If True, create a legend for the different linestyles defining
                different mass ratio of the initial concentration tracers
            legendLinestylesLoc : str
                Location of the legend for the different linestyles
        """
        assert type(ncol) is int and 0 < ncol
        assert type(simulationIds) is list
        assert type(subPlot) is bool
        assert type(axesResultSmall) is list and len(axesResultSmall) == 4
        assert type(subPlotModelYear) is int and 0 <= subPlotModelYear

        if subPlot:
            self.__axesResultSmall = plt.axes(axesResultSmall)

        #Plot the spin up for the simulationIds
        linestyles= []
        for (simulationId, label) in simulationIds:
            data = self._database.read_spinup_values_for_simid(simulationId)
            try:
                (linestyle, colorIndex, legend) = self._linestyle_color(label)
                if self._linestyleLabel[linestyle] not in linestyles:
                    linestyles.append(self._linestyleLabel[linestyle])

                if legend:
                    self._axesResult.plot(data[:,0], data[:,1], color=self._colors[colorIndex], linestyle=self._linestyle[linestyle], label=label)
                else:
                    self._axesResult.plot(data[:,0], data[:,1], color=self._colors[colorIndex], linestyle=self._linestyle[linestyle])
                if subPlot:
                    self.__axesResultSmall.plot(data[subPlotModelYear:,0], data[subPlotModelYear:,1], color=self._colors[colorIndex], linestyle=self._linestyle[linestyle])
            except IOError as e:
                print("Error message: " + e.args[1])
                print("Error message: Figure with was not created.")

        #Create legend entry for Constant
        self._axesResult.plot(np.NaN, np.NaN, color=self._colors[4], linestyle=self._linestyle[0], label='Constant')

        #Create legend for different linestyles
        if 'legendLinestyles' in kwargs and kwargs['legendLinestyles'] and len(linestyles) > 1:
            ax2 = self._axesResult.twinx()
            for ss, sty in enumerate(self._linestyle):
                ax2.plot(np.NaN, np.NaN, linestyle=self._linestyle[ss], label=self._linestyleLabel[ss], c='black')
                ax2.get_yaxis().set_visible(False)
                ax2.legend(loc=kwargs['legendLinestylesLoc'], title='Mass ratio', borderaxespad=0.2, labelspacing=0.2, borderpad=0.2, handlelength=1.0, handletextpad=0.5)

        #Set labels
        self._axesResult.set_xlabel(r'Model years [\si{{\Modelyear}}]')
        self._axesResult.set_ylabel(r'Norm [\si{\milli\mole\Phosphat\per\cubic\meter}]')
        self._axesResult.set_yscale('log', basey=10)

        if 'legend' in kwargs and kwargs['legend']:
            if 'labelspacing' in kwargs and 'handlelength' in kwargs and 'handletextpad' in kwargs and 'columnspacing' in kwargs and 'borderaxespad' in kwargs and 'borderpad' in kwargs:
                self._axesResult.legend(loc='best', ncol=ncol, labelspacing=kwargs['labelspacing'], handlelength=kwargs['handlelength'], handletextpad=kwargs['handletextpad'], columnspacing=kwargs['columnspacing'], borderaxespad=kwargs['borderaxespad'], borderpad=kwargs['borderpad'])


    def plot_tracer_norm_data(self, referenceSimulationId, simulationIds, norm='2', trajectory='', year=None, ncol=3, **kwargs):
        """
        Plot the norm of the tracer concentration difference
        
        Plot the development over the spin up (10000 years) of the difference
        between the 1dt solution and solutions calculated with different initial
        values in the norm for the given simulationIds.

        Parameters
        ----------
        referenceSimulationId : int
            SimulationId of the the reference solution.
        simulationIds : list [tuple]
            List of simulations using the simulationId and label defined in the
            tuples.
        norm : string, default: '2'
            Used norm
        trajectory : str, default: ''
            Use for '' the norm only at the first time point in a model year
            and use for 'trajectory' the norm over the whole trajectory
        year : int, default: None
            Use the reference solution (1dt solution) at the given year (e.g.
            the reference solution after a spin up over 10000 model years). If
            the value is None, use the same year for the reference and other
            solution.
        ncol : int, default: 3
            Number of columns for the legend
        **kwargs : dict
            Additional keyword arguments with keys:

            legend : bool
                If True, plot the legend
            handlelength : float
                The length of the legend handles, in font-size units.
            handletextpad : float
                The pad between the legend handle and text, in font-size units.
            columnspacing : float
                The spacing between columns, in font-size units.
            borderaxespad : float
                The pad between the axes and legend border, in font-size units.
        """
        assert type(referenceSimulationId) is int and 0 <= referenceSimulationId
        assert type(simulationIds) is list
        assert norm in DB_Constants.NORM
        assert trajectory in ['', 'trajectory']
        assert year is None or type(year) is int and 0 <= year
        assert type(ncol) is int and 0 < ncol

        #Get reference solution
        if year is None:
            data_solution = self._database.read_tracer_norm_values_for_simid(referenceSimulationId, norm=norm, trajectory=trajectory)
        else:
            data_solution = self._database.read_tracer_norm_value_for_simid_year(referenceSimulationId, year, norm=norm, trajectory=trajectory)

        #Plot the norm for the extra simulationIds
        for (simulationId, label) in simulationIds:
            (linestyle, colorIndex, legend) = self._linestyle_color(label)
            data = self._database.read_tracer_difference_norm_values_for_simid(simulationId, referenceSimulationId, yearB=year, norm=norm, trajectory=trajectory)
            if year is None:
                data[:,1] = data[:,1] / data_solution[:,1]
            else:
                data[:,1] = data[:,1] / data_solution
            try:
                self._axesResult.plot(data[:,0], data[:,1], color=self._colors[colorIndex], linestyle=self._linestyle[linestyle], label=label)

            except IOError as e:
                print("Error message: " + e.args[1])
                print("Error message: Figure with was not created.")

        #Set labels
        self._axesResult.set_xlabel(r'Model years [\si{{\Modelyear}}]')
        self._axesResult.set_ylabel(r'Relative error')
        self._axesResult.set_yscale('log', basey=10)
        if 'legend' in kwargs and kwargs['legend']:
            self._axesResult.legend(loc='best', ncol=ncol, handlelength=kwargs['handlelength'], handletextpad=kwargs['handletextpad'], columnspacing=kwargs['columnspacing'], borderaxespad=kwargs['borderaxespad'])


    def _colorMapping(self, concentrationTyp, differentTracer, distribution, tracerDistribution):
        """
        Map the initial concentration kind to color

        Parameters
        ----------
        concentrationTyp : {'vector', 'constant'}
            Use constant initial concentration or an initial concentration
            defined with vectors
        differentTracer : int
            Number of different tracer
        distribution : {None, 'Uniform', 'Normal', 'Lognormal', 'OneBox'}
            If concentrationTyp is 'vector', distribution used to generate the
            initial tracer concentration vectors.
        tracerDistribution : {None, 'set_mass', 'random_mass'}
            If 'set_mass', use the standard ratio of the mass between the
            differnt tracers and, if 'random_mass', use the random ratio between
            the different tracer

        Returns
        -------
        tuple
            Tuple including the colorIndex and the marker
        """
        assert concentrationTyp in Metos3d_Constants.METOS3D_MODEL_TRACER_CONCENTRATIONTYP
        assert differentTracer is None or type(differentTracer) is int and 1 <= differentTracer and differentTracer <= len(Metos3d_Constants.TRACER_MASK)
        assert distribution is None or distribution in Metos3d_Constants.METOS3D_MODEL_TRACER_DISTRIBUTION
        assert tracerDistribution is None or tracerDistribution in Metos3d_Constants.METOS3D_MODEL_TRACER_TRACERDISTRIBUTION

        legend = True

        #Set marker
        marker = self._marker[0]
        if concentrationTyp == 'constant' or tracerDistribution == 'random_mass':
            marker = self._marker[1]
            legend = False

        colorIndex = -1
        if concentrationTyp == 'constant':
            colorIndex = 4
        elif concentrationTyp == 'vector' and distribution == 'Uniform':
            colorIndex = 1
        elif concentrationTyp == 'vector' and distribution == 'Normal':
            colorIndex = 0
        elif concentrationTyp == 'vector' and distribution == 'Lognormal':
            colorIndex = 2
        elif concentrationTyp == 'vector' and distribution == 'OneBox':
            colorIndex = 3

        return (colorIndex, marker, legend)


    def plot_scatter_spinup_norm(self, model, initialValueParameter, year=10000, norm='2', trajectory='', alpha=0.75, **kwargs):
        """
        Scatter plot of the spin up against the norm

        Plot a scatter plot of the relation between the spin up norm and the
        norm of the tracer difference between the spin up calculation using
        different initial values and the spin up calculation using the standard
        constant inital value. The plot visualizes the ratio for the given
        initial value configuration using different colors.

        Parameters
        ----------
        model : str
            Name of the biogeochemical model
        initialValueParameter : list [tuple]
            List of tuple with the configurations of the initial value
            consisting of TODO (startTimestep, stepYear, tolerance, rho, eta, cpus,
            norm, checkConcentration, singleStep, singleStepYear, label)
        year : int, default: 10000
            Used model year of the spin up (for the spin up is the previous
            model year used)
        norm : string, default: '2'
            Used norm
        trajectory : str, default: ''
            Use for '' the norm only at the first time point in a model year
            and use for 'Trajectory' the norm over the whole trajectory
        alpha : float, default: 0.75
            The alpha blending value of the scatter plot, between 0
            (transparent) and 1 (opaque)
        **kwargs : dict
            Additional keyword arguments with keys:

            legendMarkers : bool
                If True, create legend with the different markers defining the
                different initial concentrations
            legendMarkersLoc : str
                Location of the legend for the different markers
            xscalelog : bool
                If True, use a logarithmic scale for the x axis
            xlabelpad : float
                Spacing in points from the axes bounding box including ticks
                and tick labels
            xrange : list [float]
                Restriction of the spin-up value to the given range. If a value
                is None, is range is open on this side
            yrange : list [float]
                Restriction of the relative error value to the given range. If
                a value is None, is range is open on this side
        """
        assert model in Metos3d_Constants.METOS3D_MODELS
        assert type(initialValueParameter) is list
        assert type(year) is int and 0 <= year
        assert norm in DB_Constants.NORM
        assert trajectory in ['', 'Trajectory']
        assert type(alpha) is float and 0.0 <= alpha and alpha <= 1.0

        markers = []

        for (concentrationTyp, differentTracer, distribution, tracerDistribution, label) in initialValueParameter:
            data = self._database.read_spinupNorm_relNorm_for_model_year(model, concentrationTyp=concentrationTyp, differentTracer=differentTracer, distribution=distribution, tracerDistribution=tracerDistribution, year=year, norm=norm, trajectory=trajectory, lhs=False)

            try:
                (colorIndex, marker, legend) = self._colorMapping(concentrationTyp, differentTracer, distribution, tracerDistribution)
                if marker not in markers:
                    markers.append(marker)

                #Restrict data
                x = data[:,2]
                y = data[:,3]
                if 'xrange' in kwargs and len(kwargs['xrange']) == 2:
                    if type(kwargs['xrange'][0]) is float:
                        xNew = x[x>kwargs['xrange'][0]]
                        y = y[x>kwargs['xrange'][0]]
                        x = xNew
                    if type(kwargs['xrange'][1]) is float:
                        xNew = x[x<kwargs['xrange'][1]]
                        y = y[x<kwargs['xrange'][1]]
                        x = xNew
                if 'yrange' in kwargs and len(kwargs['yrange']) == 2:
                    if type(kwargs['yrange'][0]) is float:
                        yNew = y[y>kwargs['yrange'][0]]
                        x = x[y>kwargs['yrange'][0]]
                        y = yNew
                    if type(kwargs['yrange'][1]) is float:
                        yNew = y[y<kwargs['yrange'][1]]
                        x = x[y<kwargs['yrange'][1]]
                        y = yNew

                if legend:
                    self._axesResult.scatter(x, y, s=4.0, marker=marker, color=self._colors[colorIndex], alpha=alpha, label = r'{}'.format(label))
                else:
                    self._axesResult.scatter(x, y, s=4.0, marker=marker, color=self._colors[colorIndex], alpha=alpha)
            except IOError as e:
                print("Error message: " + e.args[1])
                print("Error message: Figure was not created.")

        #Create legend entry for Constant
        self._axesResult.scatter(np.NaN, np.NaN, s=4.0, marker=self._marker[0], color=self._colors[4], label='Constant')

        #Create legend for different linestyles
        if 'legendMarkers' in kwargs and kwargs['legendMarkers'] and len(markers) > 1:
            ax2 = self._axesResult.twinx()
            for ss, sty in enumerate(markers):
                ax2.scatter(np.NaN, np.NaN, s=4.0, marker=self._marker[ss], label=self._linestyleLabel[ss], c='black')
                ax2.get_yaxis().set_visible(False)
                ax2.legend(loc=kwargs['legendMarkersLoc'], title='Mass ratio', borderaxespad=0.2, labelspacing=0.2, borderpad=0.2, handlelength=0.4, handletextpad=0.5)

        if 'xscalelog' in kwargs and kwargs['xscalelog']:
            self._axesResult.set_xscale('log', basex=10)
        self._axesResult.set_yscale('log', basey=10)
        self._axesResult.set_xlabel(r'Norm [\si{\milli\mole\Phosphat\per\cubic\meter}]', labelpad=kwargs['xlabelpad'])
        self._axesResult.set_ylabel(r'Relative error')

