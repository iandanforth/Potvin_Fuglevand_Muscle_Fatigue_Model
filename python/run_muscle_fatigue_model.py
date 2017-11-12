from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import muscle_fatigue
from yaml import load
from pprint import pprint
import numpy as np

class Params(object):
    def __init__(self, paramsDict):
        self.__dict__.update(paramsDict)

def main(arguments):
    config = None
    if arguments.config:
        config = load(arguments.config)
    else:
        print("No configuration file provided, using default config ...")
        with open('default_config.yaml') as fh:
            config = load(fh)

    if arguments.verbose:
        print("Parameters for run '%s':" % config["run_name"])
        pprint(config["parameters"])

    # Create isotonic data
    p = Params(config["parameters"])
    fthsamp = p.fthtime * p.samprate
    fth = np.zeros(fthsamp) + p.fthscale

    # Create Ramp Plateau Data
    # NOTE: SKIPPED - TODO come back to this

    # Calculations from the Fuglevand, Winter & Patla (1993) Model

    # NOTE: Ommitted lines 86/87 from direct port as they appear to do nothing

    # Recruitment Threshold Excitation (thr) for all neurons TODO: rename
    # This varies across the population
    thr = np.arange(p.nu, dtype='float64')
    n = np.arange(p.nu, dtype='float64')
    # this was modified from Fuglevand et al (1993) RTE(i) equation (1)
    # as that did not create the exact range of RTEs (ie. 'r') entered
    # Scaling factor of recruitment thresholds across neuron population
    b = np.log(p.r + (1 - p.mthr)) / (p.nu - 1)

    # Log thresholds
    thr *= b
    # Base thresholds
    thr = np.exp(thr)
    # Apply gain parameter to thresholds
    thr *= p.a
    # Shift all thresholds by minimum threshold
    thr -= 1 - p.mthr

    # Peak Firing Rate (frp) TODO: rename
    # This curve starts high with the first recruited motor unit and then trends down
    # modified from Fuglevand et al (1993) PFR equation (5) to remove thr(1)
    # before ratio
    firing_rate_diff = p.pfr1 - p.pfrL
    frp = p.pfr1 - ((thr - thr[0]) / (p.r - thr[0]) * firing_rate_diff)
    maxex = thr[p.nu -1] + (p.pfrL - p.minfr) / p.mir # maximum excitation
    maxact = int(round(maxex * p.res)) # max excitation * resolution
    ulr = (thr[p.nu -1] * 100) / maxex # recruitment threshold (%) of last motor unit

    # Calculation of the rested motor unit twitch properties (these will
    # change with fatigue)

    # Firing Rates for each MU with increased excitation (act)
    # Pre-calculate firing rates for each motor unit for each step of activation level
    # Step size is inversely proportional to the resolution parameter

    mufr = np.zeros((p.nu, maxact))
    # TODO: Convert to masked vector operations
    # act_range = np.arange(1, maxact + 1, dtype='float64')
    # # activation increment
    # print(maxex)
    # act_increment = 1 / p.res
    # print(act_increment)
    # act_range *= act_increment
    # threshold = 22.0
    # act_range -= 22
    # act_range *= p.mir
    # act_range += p.minfr
    for mu in range(p.nu):
        for act in range(maxact):
            acti = (act + 1) / p.res
            if acti >= thr[mu]:
                mufr[mu, act] = p.mir * (acti - thr[mu]) + p.minfr
                if mufr[mu, act] > frp[mu]:
                    mufr[mu, act] = frp[mu]
            else:
                mufr[mu, act] = 0.0

    k = np.arange(maxact, dtype='float64') # range of excitation levels

    # Twitch peak force (peak_twitch_forces)
    # this was modified from Fuglevand et al (1993) P(i) equation (13)
    # as that didn't create the exact range of twitch tensions (ie. 'rp') entered
    b = np.log(p.rp) / (p.nu - 1)
    peak_twitch_forces = np.exp((b * n))

    # Twitch contraction times
    c = np.log(p.rp) / np.log(p.rt)
    ct = np.zeros(p.nu)
    # assigns contraction times to each motor unit
    for mu in range(p.nu):
        # TODO: More vectorizing
        # one = 1 / peak_twitch_forces[mu]
        # print(one)
        # two = one ** (1 / c)
        # print(two)
        # three = np.dot(p.tL, two)
        # print(three)
        ct[mu] = p.tL * ((1 / peak_twitch_forces[mu]) ** (1 / c))

    # Normalized motor unit firing rates (nmufr) with increased excitation (act)
    nmufr = np.zeros((p.nu, maxact))
    for mu in range(p.nu):
        for act in range(maxact):
            nmufr[mu, act] = ct[mu] * (mufr[mu, act] / 1000)

    # Motor unit force, relative to full fusion (Pr) with increasing excitation
    # based on Figure 2 of Fuglevand et al (1993)
    sPr = 1 - np.exp((-2 * (0.4 ** 3)))
    relative_mu_forces = np.zeros((p.nu, maxact))
    for mu in range(p.nu):
        for act in range(maxact):
            # linear portion of curve
            # relative_mu_forces = MU force relative to rest 100% max excitation of 67
            if nmufr[mu, act] <= 0.4:
                relative_mu_forces[mu, act] = sPr * (nmufr[mu, act] / 0.4)
            else:
                relative_mu_forces[mu, act] = 1 - np.exp((-2 * (nmufr[mu, act] ** 3)))

    # Motor unit force with increased excitation
    mu_forces = np.zeros((p.nu, maxact))
    for mu in range(p.nu):
        for act in range(maxact):
            mu_forces[mu, act] = relative_mu_forces[mu, act] * peak_twitch_forces[mu]

    total_forces = np.sum(mu_forces, 0) # sum of forces across MUs for each excitation (dim 1)
    max_force = total_forces[-1]

    # Total force across all motor units when rested
    forces_now = np.zeros((p.nu, fthsamp))
    forces_now[:, 0] = peak_twitch_forces

    # Calculation of Fatigue Parameters (recovery currently set to zero in
    # this version)

    # note, if rp = 100 & fat = 180, there will be a 100 x 180 = 1800-fold difference in
    # the absolute fatigue of the highest threshold vs the lowest threshold.
    # The highest threshold MU will only achieve ~57# of its maximum (at 25 Hz), so the actual range of fatigue
    # rates is 1800 x 0.57 = 1026

    # fatigue rate for each motor unit
    fatigue_scale = np.log(p.fat) / (p.nu - 1)
    mu_fatigue_rates = np.exp(n * fatigue_scale)

    fatigues = np.zeros(p.nu)
    for mu in range(p.nu):
        fatigues[mu] = (mu_fatigue_rates[mu] * (p.FatFac / p.fat)) * peak_twitch_forces[mu]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c', '--config', help="Configuration file", type=argparse.FileType('r'))
    parser.add_argument('-v', '--verbose', help="Verbose output", action='store_true')
    args = parser.parse_args()
    sys.exit(main(args))
