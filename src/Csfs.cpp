// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Csfs.hpp"

#include "EigenTypes.hpp"

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <stdexcept>

namespace asmc {


CSFS::CSFS(std::map<double, CSFSEntry> CSFS_) : mCSFS(CSFS_),
  mArraySpectrum({}), mCSFS({}), mAscertainedCSFS({}),
  mFoldedAscertainedCSFS({}), mCompressedAscertainedEmissionTable({}),
  mArraySamplingFactors({}), mSamples(0) {
  if (mCSFS.size() > 0) mFoldedCSFS = foldCSFS(mCSFS);
}

CSFS::loadFromFile(std::string_view filename) {
    BufferedReader br = Utils.openFile(fileName);
    String line;
    TreeMap<Double, CSFSEntry> parsedCSFS = new TreeMap<Double, CSFSEntry>();
    while ((line = br.readLine()) != null) {
        String[] splitString = line.split("\\s+");
        if (splitString[0].compareToIgnoreCase("Time:") == 0) {
            // this is a new CSFS line entry. Expect to read time vector, size vector, Mu, Samples, Intervals, then 3 CSFS lines.
            // parse time line
            ArrayList<Double> timeVector = new ArrayList<Double>();
            for (int i = 1; i < splitString.length; i++) {
                double thisTime = Double.parseDouble(splitString[i]);
                timeVector.add(thisTime);
            }
            // parse size line
            line = br.readLine();
            splitString = line.split("\\s+");
            if (splitString[0].compareToIgnoreCase("Size:") != 0) {
                Utils.exit("Parsed line \"" + line + "\" after Time line in " + fileName + ".");
            }
            ArrayList<Double> sizeVector = new ArrayList<Double>();
            for (int i = 1; i < splitString.length; i++) {
                double thisSize = Double.parseDouble(splitString[i]);
                sizeVector.add(thisSize);
            }
            // parse mu
            line = br.readLine();
            splitString = line.split("\\s+");
            if (splitString[0].compareToIgnoreCase("Mu:") != 0) {
                Utils.exit("Parsed line \"" + line + "\" after Size line in " + fileName + ".");
            }
            double mu = Double.parseDouble(splitString[1]);
            // parse samples
            line = br.readLine();
            splitString = line.split("\\s+");
            if (splitString[0].compareToIgnoreCase("Samples:") != 0) {
                Utils.exit("Parsed line \"" + line + "\" after Mu line in " + fileName + ".");
            }
            int samples = Integer.parseInt(splitString[1]);
            // parse samples
            line = br.readLine();
            splitString = line.split("\\s+");
            if (splitString[0].compareToIgnoreCase("Interval:") != 0) {
                Utils.exit("Parsed line \"" + line + "\" after Samples line in " + fileName + ".");
            }
            double from = Double.parseDouble(splitString[1]);
            double to = Double.parseDouble(splitString[2]);
            // parse CSFS
            double[][] CSFS = new double[3][samples - 1];
            for (int dist = 0; dist < 3; dist++) {
                line = br.readLine();
                splitString = line.split("\\s+");
                for (int undist = 0; undist < splitString.length; undist++) {
                    double thisEntry = Double.parseDouble(splitString[undist]);
                    CSFS[dist][undist] = thisEntry;
                }
            }
            CSFSEntry thisEntry = new CSFSEntry(timeVector, sizeVector, mu, from, to, samples, CSFS);
            parsedCSFS.put(from, thisEntry);
        } else {
            Utils.exit("Badly formatted CSFS file. Parsed line \"" + line + "\", which does not contain a Time definition. ");
        }
    }
    Utils.print("Read " + parsedCSFS.size() + " CSFS entries.");
    return new CSFS(parsedCSFS);
}

public boolean verify(ArrayList<Double> timeVectorOriginal, ArrayList<Double> sizeVectorOriginal, double mu, int samples, ArrayList<Double> discretizationOriginal) {
    try {
        ArrayList<Double> timeVector = (ArrayList<Double>) timeVectorOriginal.clone();
        timeVector.remove(timeVector.size() - 1);
        ArrayList<Double> sizeVector = (ArrayList<Double>) sizeVectorOriginal.clone();
        sizeVector.remove(sizeVector.size() - 1);
        ArrayList<Double> discretization = (ArrayList<Double>) discretizationOriginal.clone();
        discretization.remove(discretization.size() - 1);
        for (double from : discretization) {
            if (!CSFS.containsKey(from)) {
                Utils.warning("CSFS does not contain interval " + from + ".");
                return false;
            }
            CSFSEntry thisEntry = CSFS.get(from);
            if (thisEntry.mu != mu) {
                Utils.warning("CSFS entry " + from + " has different mu: " + thisEntry.mu + ".");
                return false;
            }
            if (!compareArrays(thisEntry.timeVector, timeVector)) {
                Utils.print(thisEntry.timeVector);
                Utils.print(timeVector);
                Utils.warning("CSFS entry " + from + " has different time vector.");
                return false;
            }
            if (!compareArrays(thisEntry.sizeVector, sizeVector)) {
                Utils.warning("CSFS entry " + from + " has different size vector.");
                return false;
            }
            if (thisEntry.samples != samples) {
                if (samples == Integer.MAX_VALUE) {
                    samples = thisEntry.samples;
                } else {
                    Utils.warning("CSFS entry " + from + " has different samples (want: " + samples + ", found: " + thisEntry.samples + ")");
                    return false;
                }
            }
        }
        return true;
    } catch (Exception e) {
        Utils.warning("Something went wrong.");
        return false;
    }

}

@Override
public String toString() {
    StringBuilder sb = new StringBuilder();
    for (double from : CSFS.keySet()) {
        sb.append(CSFS.get(from));
    }
    return sb.toString();
}

void CSFS::fixAscertainment(Data data, int samples, Transition transition) {
    computeArraySamplingFactors(data, samples, transition);
    // CSFS is loaded here, but fixed later.
    for (std::pair<double, CSFSEntry> entry: mCSFS) mAscertainedCSFS.emplace(entry);
    applyFactors();
    mFoldedAscertainedCSFS = foldCSFS(mAscertainedCSFS);
    mCompressedAscertainedEmissionTable = compressCSFS(mFoldedAscertainedCSFS);
}

mat_dt CSFS::computeClassicEmission(array_dt expectedTimes, double mu) {
  vec_dt emissionRow;
  mat_dt emission;
  for(double interval : mExpectedTimes) emissionRow << std::exp(-2 * interval * mu);
  emission.resize(2, emissionRow.size());
  emission.row(0) = emissionRow;
  emission.row(1) = 1 - emissionRow;
  return emission;
}

void CSFS::computeArraySamplingFactors(Data data, int samples, Transition transition) {
    mSamples = samples;
    std::vector<double> coalDist = transition.getCoalDist();
    std::vector<double> AFS;
    // double[] AFS = new double[samples];
    // the first entry of the CSFS may not be zero, since it's a shared doubleton
    int counter = 0;
    for (double from : CSFS.keySet()) {
        double coalescentProbabilityThisTimeSlice = coalDist[counter];
        CSFSEntry thisCsfsDM = CSFS.get(from);
        for (int row = 0; row < 3; row++) {
            for (int column = 0; column < samples - 1; column++) {
                int pos = row + column;
                if (pos > samples / 2) {
                    pos = samples - pos;
                }
                AFS[pos] += coalescentProbabilityThisTimeSlice * thisCsfsDM.CSFS[row][column];
            }
        }
        counter++;
    }

    // normalize spectrum
    AFS[0] = 0.;
    double norm = 0.;
    for (int i = 0; i < samples; i++) {
        norm += AFS[i];
    }
    for (int i = 0; i < samples; i++) {
        AFS[i] /= norm;
    }

    // fold AFS
    int halfTotal = samples / 2;
    for (int i = halfTotal + 1; i < samples; i++) {
        AFS[samples - i] += AFS[i];
        AFS[i] = 0;
    }
    // normalize spectrum
    norm = 0.;
    for (int i = 0; i < samples; i++) {
        norm += AFS[i];
    }
    for (int i = 0; i < samples; i++) {
        AFS[i] /= norm;
    }

    // foldedAFS contains probability a site has MAF i given the site is polymorphic in the sequence data
    double[] foldedAFS = new double[halfTotal + 1];
    for (int i = 0; i <= halfTotal; i++) {
        foldedAFS[i] = AFS[i];
    }

    // now get foldedAFS_array, the probability a site has MAF i given it is polymorphic in the sample (array)
    this.arraySpectrum = new ArraySpectrum(data, samples);
    double[] foldedAFS_array = arraySpectrum.spectrum;
    for (int i = 0; i < foldedAFS_array.length; i++) {
    }
    double[] samplingFactors = new double[halfTotal + 1];

    for (int i = 1; i < foldedAFS_array.length; i++) {
        samplingFactors[i] = foldedAFS_array[i] / foldedAFS[i];
    }
    arraySamplingFactors = samplingFactors;
}

private void applyFactors() {
    // apply sampling factors and renormalize
    // note that the first entry of the CSFS may not be zero, since it's a shared doubleton
    for (double from : ascertainedCSFS.keySet()) {
        double[][] thisCSFS = ascertainedCSFS.get(from).CSFS;
        thisCSFS[0][0] = 0.;
        double norm = 0.;
        for (int row = 0; row < 3; row++) {
            for (int column = 0; column < samples - 1; column++) {
                // if the spectrum is folded, this emission is mapped to this position
                int pos = row + column;
                if (pos > samples / 2) {
                    pos = samples - pos;
                }
                // and if we're looking at array data, this MAF is adjusted using this factor
                double factor = arraySamplingFactors[pos];
                // apply factor
                thisCSFS[row][column] *= factor;
                // sum value to renomralize to 1 later on
                norm += thisCSFS[row][column];
            }
        }
        norm /= 1 - arraySpectrum.monomorphic;
        for (int row = 0; row < 3; row++) {
            for (int column = 0; column < samples - 1; column++) {
                thisCSFS[row][column] /= norm;
            }
        }
        thisCSFS[0][0] = arraySpectrum.monomorphic;
        ascertainedCSFS.get(from).CSFS = thisCSFS;
    }
}

public TreeMap<Double, CSFSEntry> foldCSFS(TreeMap<Double, CSFSEntry> CSFS) {
    TreeMap<Double, CSFSEntry> foldedCSFS = new TreeMap<Double, CSFSEntry>();
    int samples = CSFS.firstEntry().getValue().samples;
    int undistinguished = samples - 2;
    for (double from : CSFS.keySet()) {
        CSFSEntry foldedEntry = new CSFSEntry(CSFS.get(from));
        double[][] thisCsfs_double = foldedEntry.CSFS;
        // code to fold the spectrum
        if (samples % 2 != 0) {
            Utils.exit("ConditionalSFS called with odd number of samples.");
        }
        int half = samples / 2;
        double[][] thisCsfs_double_folded = new double[2][half + 1];
        for (int row = 0; row < 3; row++) {
            for (int column = 0; column < undistinguished + 1; column++) {
                Integer[] coord = new Integer[]{row, column};
                Integer[] foldedCoord = getFoldedObservationFromUnfolded(coord, samples);
                thisCsfs_double_folded[foldedCoord[0]][foldedCoord[1]] += thisCsfs_double[row][column];
            }
        }
        foldedEntry.CSFS = thisCsfs_double_folded;
        foldedCSFS.put(from, foldedEntry);
    }
    return foldedCSFS;
}

std::pair<int, int> CSFS::getFoldedObservationFromUnfolded(std::pair<int, int> unfolded, int totalSamples) {
  auto [dist, undist] = unfolded;
  if (totalSamples % 2 != 0) throw std::runtime_error(
      "Function getFoldedObservationFromUnfolded was called with odd total sample size. "
      "Only diploid samples are supported at the moment.");
  if (undist + dist > totalSamples / 2) undist = totalSamples - 2 - undist; // flip
  if (dist == 2) dist = 0;
  return std::make_pair(dist, undist);
}

public double[][] compressCSFS(TreeMap<Double, CSFSEntry> CSFS) {
    double[][] compressed = new double[2][CSFS.size()];
    int timeInterval = 0;
    for (double from : CSFS.keySet()) {
        double[][] thisCSFS = CSFS.get(from).CSFS;
        for (int k = 0; k < thisCSFS[0].length; k++) {
            compressed[0][timeInterval] += thisCSFS[0][k];
            compressed[1][timeInterval] += thisCSFS[1][k];
        }
        timeInterval++;
    }
    return compressed;
}

} // namespace asmc
