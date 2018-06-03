#! /usr/bin/env python

from __future__ import division
# import modules needed
import math
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
# Files imported
from Cylinder_Geometry import Cylinder_Geometry
# import SimFunctions #not used anymore



############
# References#
############
# 1. Internal Combustion Engines: Applied Thermal Sciences 2nd edition
# author(s): Ferguson, Colin R.; Kirkpatrick, Allan T.
# 2.

# NOTES!
# All units are SI standard (ex. kelvin, kilogram, meter, second, joules, pascals)

class System:
    def __init__(self):
        # properties
        self.ambient_temp = 300  # K
        self.ambient_pressure = 101325  # Pa

        # self.composition
        self.geom = Cylinder_Geometry()
        self.bore = 0.095  # m
        self.stroke = 0.0634  # m
        self.connecting_rod_length = 0.1074  # m
        self.compression_ratio = 13.5
        self.swept_volume = 0.0004494  # m^3
        self.clearance_volume = 0.00003808  # m^3
        self.total_volume = self.swept_volume + self.clearance_volume

        # intake
        self.intake_diameter = 0.02
        self.intake_area = math.pi * (self.intake_diameter ** 2) / 4

        # Valve Timing
        self.rpm = 5000
        self.intake_open = 40
        self.intake_close = 120
        self.exhaust_open = 535
        self.exhaust_close = 0
        self.spark_advance = -55
        self.ignition_start = 360 + self.spark_advance
        self.combustion_duration = 40
        self.AFR = 14.5  # mass ratio
        self.volumetric_efficiency = 0.8

        # fuel
        self.fuelType = 'gasoline'
        self.LHVgasoline = 44.40e6  # J/kg
        self.air_density = 1.225  # kg/m^3
        #	#q_combustion appears off by 1000 times...units?
        self.Q_combustion = self.volumetric_efficiency * self.swept_volume * self.air_density * (
                1 / self.AFR) * self.LHVgasoline

        # setup lists
        self.crank_angle = []
        self.volume = []
        self.time = []
        self.pressure = [self.ambient_pressure]  ### THIS NEEDS TO HAVE AN INTAKE SIM
        self.temp = [self.ambient_temp]
        self.intake_valve = [0]
        self.exhaust_valve = [0]
        self.xb = [0]

        # discretization
        self.step = 0.1

    #############functions#####################
    def instantVolume(self, crank_angle):
        # cylinder volume at specific crank angle
        #print (crank_angle)
        instVol = self.clearance_volume + ((1 - self.geom.piston_position(crank_angle)) * self.swept_volume)
        return instVol

    def angle2time(self, angle_1, angle_2):
        # time difference between two crank angles
        dt = (angle_2 - angle_1) / (self.rpm * 60)
        return dt

    def gammaAir(self, temp):
        # gamma of air at specific temperature
        # source...
        return 1.40 - (7.18e-5 * temp)

    def airProperties(self, property):
        # properties of air at ~278 K
        self.basic_properties = {
            'gamma': 1.4,
            'cv': 0.714,
            'cp': 1,
            'R': 287.058,
            'density': 1.225
        }
        return self.basic_properties[property]

    def areaCylinder(self, crank_angle):
        # surface area of cylinder at instant
        piston_area = math.pi * (self.bore ** 2) * 0.25
        wall_area = math.pi * self.bore * self.stroke * (1 - self.geom.piston_position(crank_angle))
        return wall_area + (piston_area * 2)

    # heat transfer coefficient function needs checking out, its throwing errors
    def heatTransferCoefficientWoschni(self, dQin=0):
        # equation 8.39 from Ferguson
        mean_piston_speed = 2 * self.stroke * self.rpm * (1 / 60)
        if self.intake_valve[-1] > 0 or self.exhaust_valve[-1] > 0:
            U = 6.18 * mean_piston_speed
        elif self.intake_valve[-1] == 0 and self.exhaust_valve[-1] == 0:
            V1 = self.volume[-2]
            gam = self.gammaAir(self.temp[-1])
            # state properties at intake valve closing
            P0 = self.state_intake_close[0]
            T0 = self.state_intake_close[1]
            V0 = self.state_intake_close[2]
            Vd = self.swept_volume  # displacement volume
            dPc = ((gam - 1) / V1) * dQin  # instantaneous pressure rise due to combustion
            U = (2.28 * mean_piston_speed) + (0.00324 * T0 * (Vd * dPc) * (1 / (V0 * P0)))
        P = self.pressure[-1] * 0.001  # convert to KPa for this instance
        T = self.temp[-1]
        b = self.bore
        # convective heat transfer coefficient (woschni)
        hg = 3.26 * (P ** 0.8) * (U ** 0.8) * (b ** -0.2) * (T ** -0.55)
        return hg

    def woschniHeatLoss(self, dQin=0):
        # equation 8.27 from Ferguson
        d_theta = self.crank_angle[-1] - self.crank_angle[-2]
        # convective heat transfer coefficient (function of theta)
        # hg=self.heatTransferCoefficientWoschni(dQin)
        hg = 3000  # fix function so it's not hard coded
        Aw = self.areaCylinder(self.crank_angle[-1])
        Tgas = self.temp[-1]
        Twall = 380  # 380 kelvin=225 Fahrenheit
        N = self.rpm / 120  # convert to power cycles per second
        # delta Q transferred to cylinder wall through convection
        dQout = d_theta * hg * Aw * (Tgas - Twall) / N
        return dQout

    def pressureUpdate(self, dQin=0, dQw=0):
        # equation 8.26 from Ferguson
        V2 = self.volume[-1]
        V1 = self.volume[-2]
        dV = V2 - V1
        P1 = self.pressure[-1]
        T1 = self.temp[-1]
        gam = self.gammaAir(T1)
        return P1 + (((gam - 1) / V1) * dQin) - dQw - (gam * (P1 / V1) * (dV))

    def tempUpdate(self, dQin=0, dQw=0):
        # Ferguson
        # 2 equations together
        V2 = self.volume[-1]
        V1 = self.volume[-2]
        T1 = self.temp[-1]
        Cv = 0.714  # calculate from JANAF or NASA coefficients (buttsworth matlab files)
        gam = self.gammaAir(T1)
        return (T1 * ((V1 / V2) ** (gam - 1))) + ((dQin - dQw) / Cv)

    def combustionStep(self, theta):
        # equation 2.22 from Ferguson (Single Wiebe)
        a = 6
        m = 1.75  # these constants can change alot
        combustion_duration = self.combustion_duration
        combustion_start = self.ignition_start
        xb = 1 - math.exp(-a * ((((theta) - combustion_start) / combustion_duration) ** (m + 1)))
        return xb

    #############Simulation#####################
    def simulation(self, rpm):
        # self explanatory

        self.sim_start_time = timer()

        self.rpm = rpm
        step = self.step
        theta = np.arange(0, 720, step)

        for i in theta:
            # Calculated for each step
            self.crank_angle.append(i)
            self.volume.append(self.instantVolume(i))
            self.time.append(self.angle2time(0, i))

            #########
            # intake stroke
            # equations (keeping it simple for now)
            if i > 0 and i <= 180:
                self.pressure.append(self.ambient_pressure)
                self.temp.append(self.ambient_temp)
                self.intake_valve.append(1)
                self.exhaust_valve.append(0)

                if i == self.intake_close:
                    # saves state properties at intake valve closing for characteristic...
                    # gas velocity used in woschni heat transfer coefficient calc above.
                    self.state_intake_close = [self.pressure[-1], self.temp[-1], self.volume[-1]]

            ##########
            # compression stroke
            elif i > 180 and i <= 360:
                # valves still closed
                self.intake_valve.append(0)
                self.exhaust_valve.append(0)

                # compression w/o heat addition
                if i < self.ignition_start:
                    # update state properties
                    dQin = 0
                    self.temp.append(self.tempUpdate(dQin))
                    self.pressure.append(self.pressureUpdate(dQin))

                    if i == self.intake_close:
                        # saves state properties at intake valve closing for characteristic...
                        # gas velocity used in woschni heat transfer coefficient calc above.
                        self.state_intake_close = [self.pressure[-1], self.temp[-1], self.volume[-1]]

                # compression w/heat addition
                elif i >= self.ignition_start:
                    # calculate burnt fraction
                    self.xb.append(self.combustionStep(i))
                    dxb = self.xb[-1] - self.xb[-2]
                    dQin = dxb * self.Q_combustion  # delta heat added to system through combustion
                    # update state properties
                    self.temp.append(self.tempUpdate(dQin))
                    self.pressure.append(self.pressureUpdate(dQin))

                else:
                    print ('Mistake in compression')

            ###############
            # Expansion Stroke
            elif i > 360 and i <= 540:
                # Convective heat transfer to walls
                dQw = self.woschniHeatLoss()

                # expansion w/heat addition
                if i > 360 and i <= (self.ignition_start + self.combustion_duration):
                    # calculate burnt fraction
                    self.xb.append(self.combustionStep(i))
                    dxb = self.xb[-1] - self.xb[-2]
                    dQin = dxb * self.Q_combustion  # delta heat added to system through combustion
                    # update state properties
                    self.temp.append(self.tempUpdate(dQin, dQw))
                    self.pressure.append(self.pressureUpdate(dQin))
                    self.intake_valve.append(0)
                    self.exhaust_valve.append(0)

                # expansion w/o heat addition
                elif i > (self.ignition_start + self.combustion_duration) and i <= self.exhaust_open:
                    dQin = 0
                    # update state properties
                    self.temp.append(self.tempUpdate(dQin, dQw))
                    self.pressure.append(self.pressureUpdate(dQin))
                    self.intake_valve.append(0)
                    self.exhaust_valve.append(0)

                # expansion w/ exhaust valve open
                elif i > self.exhaust_open and i <= 540:
                    dQin = 0
                    # update state properties
                    self.temp.append(self.tempUpdate(dQin, dQw))
                    self.pressure.append(self.pressureUpdate(dQin))
                    self.intake_valve.append(0)
                    self.exhaust_valve.append(1)

                else:
                    print ('mistake in expansion')

            ###############
            # Exhaust Stroke
            elif i > 540 and i <= 720:
                self.temp.append(self.temp[-1])
                self.pressure.append(self.pressure[-1])
                self.intake_valve.append(0)
                self.exhaust_valve.append(1)

        self.sim_end_time = timer()
        return self.brakeTorque() * -0.1

    ######################################
    # simulation data processing functions

    def plotPV(self):
        plt.plot(self.volume, self.pressure)
        plt.title('P-V Diagram')
        plt.xlabel('Volume [m^3]')
        plt.ylabel('Pressure [Pa]')
        plt.show()

    def plotIt(self, yoMama):
        # temp, pres, volume, crank angle
        if yoMama == 'pressure':
            plt.plot(self.crank_angle, self.pressure)
            plt.plot([360, 360], [0, max(self.pressure) + 1E6])
            plt.title('pressure')
            plt.xlabel('crank angle [degrees]')
            plt.ylabel('pressure [Pa]')
        elif yoMama == 'temperature':
            plt.plot(self.crank_angle, self.temp)
            plt.title('temperature')
            plt.xlabel('crank angle [degrees]')
            plt.ylabel('Temperature [K]')
        elif yoMama == 'volume':
            plt.plot(self.crank_angle, self.volume)
            plt.title('volume')
            plt.xlabel('crank angle [degrees]')
            plt.ylabel('Volume [m^3]')
        elif yoMama == 'intake_valve':
            plt.plot(self.crank_angle, self.intake_valve)
        else:
            return "Sumthin Wrong"
        plt.show()

    def timeElapsed(self):
        return self.sim_end_time - self.sim_start_time

    def IMEPnet(self):
        return np.mean(self.pressure)

    def FMEP(self):
        return 200000

    def BMEP(self):
        return self.IMEPnet() - self.FMEP()

    def brakeTorque(self):
        # units [N*m]=[J]
        return self.BMEP() * self.swept_volume

    def power(self, units='watts'):
        # units [Watts]=[J/s]=[Nm/s]
        power = self.brakeTorque() * self.rpm * (1 / 60) * 2 * math.pi
        if units == 'W' or units == 'w' or units == 'Watts' or units == 'watts':
            return power
        elif units == 'KW' or units == 'kw' or units == 'Kw' or units == 'kW' or units == 'Kilowatts' or units == 'kilowatts':
            return power / 1000
        else:
            print ('Specify units, use watts or kilowatts.')


class RunSim:
    def __init__(self):
        cv = System()

    def runIt(self, range):
        l = np.arange(range[0], range[1], range[2])

        self.p = []
        self.t = []
        for i in l:
            self.t.append(cv.simulation(i))

        # plt.plot(l,self.p)
        plt.plot(l, self.t)
        plt.show()


if __name__ == "__main__":
    print('nothing yet')
    i=9000
    #for i in range(1000,12000,1000):
    geom = Cylinder_Geometry
    # print geom.piston_position(0)
    cv = System()
    # print cv.gammaAir(300)
    cv.simulation(i)
    # print len(cv.crank_angle)
    # print len(cv.volume)
    # print len(cv.temp)
    # print len(cv.pressure)
    # print cv.timeElapsed()
    cv.plotPV()
    ##cv.plotIt('temperature')
    cv.plotIt('pressure')
    # cv.plotIt('volume')
    print (cv.IMEPnet())
    # print cv.BMEP()
    print ("Brake Torque at", i, "rmp is", cv.brakeTorque(), "N-m")
    print ("Power at", i, "rmp is", cv.power(), "watts")

# rs=RunSim()
# range=[2000,8000,100]
# print rs.runIt(range)

# things to add
# calculator to find cv at each temp (see buttsworth.py and pdf)
# JANAF, NASA
# add heat loss to wall (woschni or annand)
# friction losses
# volumetric efficiency
# ignition duration/combustion duration(rate time constant)
# use Arrhenius equation
