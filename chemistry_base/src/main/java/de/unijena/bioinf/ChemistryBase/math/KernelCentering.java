package de.unijena.bioinf.ChemistryBase.math;

public class KernelCentering {

    protected final double[] averages;
    protected final double[] diagonal;
    protected final double average;

    public KernelCentering(double average, double[] averages, double[] diagonalAfterCentering) {
        this.average = average;
        this.averages = averages;
        this.diagonal = diagonalAfterCentering;
    }

    public KernelCentering(double[][] quadraticMatrix, boolean standardize) {
        final double[] rowMeans = new double[quadraticMatrix.length];
        final double[] colMeans = new double[quadraticMatrix.length];
        double mean = 0;
        for (int row=0; row < quadraticMatrix.length; ++row) {
            double m=0d;
            final double[] krow = quadraticMatrix[row];
            for (int i=0; i < krow.length; ++i) {
                m += krow[i];
                colMeans[i] += krow[i];
            }
            mean += m;
            m /= krow.length;
            rowMeans[row] = m;
        }
        for (int k=0; k < colMeans.length; ++k) {
            colMeans[k] /= quadraticMatrix.length;
        }

        mean /= (quadraticMatrix.length*quadraticMatrix[0].length);

        this.average = mean;
        this.averages = colMeans;

        if (standardize) {
            diagonal = new double[quadraticMatrix.length];
            for (int i=0; i < quadraticMatrix.length; ++i) {
                diagonal[i] = quadraticMatrix[i][i] - 2*averages[i] + average;
            }
        } else {
            diagonal = null;
        }
    }

    public KernelCentering withoutNormalizing() {
        return new KernelCentering(average,averages,null);
    }

    public double[] getMainDiagonalAfterCentering() {
        return diagonal;
    }

    public double[] getColumnAverages() {
        return averages;
    }

    public double getMatrixAverage() {
        return average;
    }

    public void applyToTrainMatrix(double[][] kernel) {
        for (int i=0; i < kernel.length; ++i) {
            for (int j=0; j <= i; ++j) {
                kernel[i][j] = kernel[i][j] - averages[i] - averages[j] + average;
                if (diagonal!=null) kernel[i][j] = MatrixUtils.norm(kernel[i][j], diagonal[i], diagonal[j]);
                if (Double.isNaN(kernel[i][j])) {
                    System.err.println("KC: " + kernel[j][i] + " becomes NaN! diagonals: " + diagonal[i] + " and " + diagonal[j]);
                }
                kernel[j][i] = kernel[i][j];
            }
        }
    }

    public void applyToKernelMatrix(double[][] kernel, double[] norms) {
        for (int row=0; row < kernel.length; ++row) {
            applyToKernelRow(kernel[row], norms==null ? Double.NaN : norms[row]);
        }
    }

    public void applyToKernelRow(double[] kernel, double norm) {
        double avg=0d;
        for (int col=0; col < kernel.length; ++col) {
            avg += kernel[col];
        }
        avg /= kernel.length;
        for (int col=0; col < kernel.length; ++col) {
            kernel[col] = kernel[col] - avg - averages[col] + average;
            if (diagonal!=null) {
                kernel[col] = MatrixUtils.norm(kernel[col], diagonal[col], norm -2 * avg + average);
            }
        }
    }

}
