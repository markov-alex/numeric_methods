package ru.mai.chm.lab2;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Supplier;

public class Solver {

    private static final Gson gson = new GsonBuilder()
            .setPrettyPrinting()
            .create();

    private final double leftBound;
    private final double rightBound;
    private final int N;
    private final int K;
    private final double T;
    private final double alpha;
    private final double beta;
    private final double gamma;
    private final double delta;
    private final BiFunction<Double, Double, Double> analyticalSolution;
    private final Function<Double, Double> phiLeft;
    private final Function<Double, Double> phiRight;
    private final Supplier<Double> psi1;
    private final Function<Double, Double> psi2;
    private final Supplier<Double> psi1SecondDerivative;

    public Solver(double leftBound, double rightBound, int n, int k, double t, double alpha,
                  double beta, double gamma, double delta, BiFunction<Double, Double, Double> analyticalSolution,
                  Function<Double, Double> phiLeft, Function<Double, Double> phiRight, Supplier<Double> psi1,
                  Function<Double, Double> psi2, Supplier<Double> psi1SecondDerivative) {
        this.leftBound = leftBound;
        this.rightBound = rightBound;
        N = n;
        K = k;
        T = t;
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
        this.delta = delta;
        this.analyticalSolution = analyticalSolution;
        this.phiLeft = phiLeft;
        this.phiRight = phiRight;
        this.psi1 = psi1;
        this.psi2 = psi2;
        this.psi1SecondDerivative = psi1SecondDerivative;
    }

    public String answerAnalyticalSolution() {
        double tau = T / K;
        double h = (rightBound - leftBound) / N;

        List<List<Double>> u = new ArrayList<>(K + 1);
        for (int i = 0; i <= K; ++i) {
            u.add(new ArrayList<>(Collections.nCopies(N + 1, 0.)));
        }

        for (int k = 0; k <= K; ++k) {
            double t = k * tau;
            for (int j = 0; j <= N; ++j) {
                double x = leftBound + j * h;
                u.get(k).set(j, analyticalSolution.apply(x, t));
            }
        }

        Answer answer = new Answer(leftBound, rightBound, tau, h, T, u);
        return gson.toJson(answer);
    }

    public String explicitScheme(ApproximationInitialConditionType approximationInitialConditionType,
                                 ApproximationBoundaryConditionType approximationBoundaryConditionType) {
        double tau = T / K;
        double h = (rightBound - leftBound) / N;

        List<List<Double>> u = new ArrayList<>(K + 1);
        for (int i = 0; i <= K; ++i) {
            u.add(new ArrayList<>(Collections.nCopies(N + 1, 0.)));
        }

        // u(x, t0) - нулевой временной слой
        for (int j = 0; j <= N; ++j) {
            u.get(0).set(j, psi1.get());
        }
        // u(x, t1) - первый временной слой
        for (int j = 0; j <= N; ++j) {
            double x = leftBound + j * h;
            if (approximationInitialConditionType == ApproximationInitialConditionType.FIRST_DEGREE) {
                u.get(1).set(j, tau * psi2.apply(x) + psi1.get());
            } else {
                u.get(1).set(j, psi1.get() * (1. - 3. * Math.pow(tau, 2) / 2.) +
                        psi2.apply(x) * tau + psi1SecondDerivative.get() * Math.pow(tau, 2) / 2.);
            }
        }

        for (int k = 1; k < K; ++k) {
            double t = (k + 1) * tau;

            for (int j = 1; j < N; ++j) {
                u.get(k + 1).set(j, u.get(k).get(j + 1) * Math.pow(tau, 2) / Math.pow(h, 2) +
                        u.get(k).get(j) * (-2. * Math.pow(tau, 2) / Math.pow(h, 2) - 3. * Math.pow(tau, 2) + 2.) +
                        u.get(k).get(j - 1) * Math.pow(tau, 2) / Math.pow(h, 2) - u.get(k - 1).get(j));
            }

            if (approximationBoundaryConditionType == ApproximationBoundaryConditionType.TWO_POINT_FIRST_DEGREE) {
                u.get(k + 1).set(0, phiLeft.apply(t) / (beta - alpha / h) - alpha / h *
                        u.get(k + 1).get(1) / (beta - alpha / h));
                u.get(k + 1).set(N, phiRight.apply(t) / (gamma / h + delta) + gamma / h *
                        u.get(k + 1).get(N - 1) / (gamma / h + delta));
            } else if (approximationBoundaryConditionType == ApproximationBoundaryConditionType.TWO_POINT_SECOND_DEGREE) {
                u.get(k + 1).set(0, (phiLeft.apply(t) - u.get(k + 1).get(1) * alpha / h +
                        (u.get(k - 1).get(0) - 2. * u.get(k).get(0)) * h * alpha / 2. / Math.pow(tau, 2)) /
                        (beta - alpha * (1. / h + 3. * h / 2. + h / 2. / Math.pow(tau, 2))));
                u.get(k + 1).set(N, (phiRight.apply(t) + u.get(k + 1).get(N - 1) * gamma / h +
                        (2. * u.get(k).get(N) - u.get(k - 1).get(N)) * h * gamma / 2. / Math.pow(tau, 2)) /
                        (delta + gamma * (1. / h + 3. * h / 2. + h / 2. / Math.pow(tau, 2))));
            } else {
                double tmp = beta - 3. * alpha / 2. / h;
                u.get(k + 1).set(0, phiLeft.apply(t) / tmp - 2. * alpha / h / tmp * u.get(k + 1).get(1) +
                        alpha / 2. / h / tmp * u.get(k + 1).get(2));
                tmp = 3. * gamma / 2. / h + delta;
                u.get(k + 1).set(N, phiRight.apply(t) / tmp + 2. * gamma / h / tmp + u.get(k + 1).get(N - 1) -
                        gamma / 2. / h / tmp * u.get(k + 1).get(N - 2));
            }
        }

        Answer answer = new Answer(leftBound, rightBound, tau, h, T, u);
        return gson.toJson(answer);
    }

    public List<Double> tridiagonalAlgo(List<Double> a, List<Double> b, List<Double> c, List<Double> d) {
        int n = a.size();
        List<Double> P = new ArrayList<>(Collections.nCopies(n + 1, 0.));
        List<Double> Q = new ArrayList<>(Collections.nCopies(n + 1, 0.));
        for (int i = 0; i < n; ++i) {
            P.set(i + 1, -1 * c.get(i) / (b.get(i) + a.get(i) * P.get(i)));
            Q.set(i + 1, (d.get(i) - a.get(i) * Q.get(i)) / (b.get(i) + a.get(i) * P.get(i)));
        }
        List<Double> x = new ArrayList<>(Collections.nCopies(n, 0.));
        for (int i = n - 1; i >= 0; --i) {
            if (i == n - 1) {
                x.set(i, Q.get(n));
            } else {
                x.set(i, Q.get(i + 1) + P.get(i + 1) * x.get(i + 1));
            }
        }
        return x;
    }

    public String implicitScheme(ApproximationInitialConditionType approximationInitialConditionType,
                                 ApproximationBoundaryConditionType approximationBoundaryConditionType) {
        double tau = T / K;
        double h = (rightBound - leftBound) / N;

        List<List<Double>> u = new ArrayList<>(K + 1);
        for (int i = 0; i <= K; ++i) {
            u.add(new ArrayList<>(Collections.nCopies(N + 1, 0.)));
        }

        // u(x, t0) - нулевой временной слой
        for (int j = 0; j <= N; ++j) {
            u.get(0).set(j, psi1.get());
        }
        // u(x, t1) - первый временной слой
        for (int j = 0; j <= N; ++j) {
            double x = leftBound + j * h;
            if (approximationInitialConditionType == ApproximationInitialConditionType.FIRST_DEGREE) {
                u.get(1).set(j, tau * psi2.apply(x) + psi1.get());
            } else {
                u.get(1).set(j, psi1.get() * (1. - 3. * Math.pow(tau, 2) / 2.) +
                        psi2.apply(x) * tau + psi1SecondDerivative.get() * Math.pow(tau, 2) / 2.);
            }
        }

        for (int k = 1; k < K; ++k) {
            double sigma = Math.pow(tau, 2) / Math.pow(h, 2);
            double t = (k + 1) * tau;
            List<Double> aCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));
            List<Double> bCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));
            List<Double> cCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));
            List<Double> dCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));

            for (int j = 1; j <= N - 1; ++j) {
                aCoefficients.set(j, 1.);
                bCoefficients.set(j, -2. - 1. / sigma - 3. * Math.pow(h, 2));
                cCoefficients.set(j, 1.);
                dCoefficients.set(j, (u.get(k - 1).get(j) - 2. * u.get(k).get(j)) / sigma);
            }

            if (approximationBoundaryConditionType == ApproximationBoundaryConditionType.TWO_POINT_FIRST_DEGREE) {
                bCoefficients.set(0, beta - alpha / h);
                cCoefficients.set(0, alpha / h);
                dCoefficients.set(0, phiLeft.apply(t));
                aCoefficients.set(N, -1. * gamma / h);
                bCoefficients.set(N, gamma / h + delta);
                dCoefficients.set(N, phiRight.apply(t));
            } else if (approximationBoundaryConditionType == ApproximationBoundaryConditionType.TWO_POINT_SECOND_DEGREE) {
                bCoefficients.set(0, beta - alpha * (1. / h + 3. * h / 2. + h / 2. / Math.pow(tau, 2)));
                cCoefficients.set(0, alpha / h);
                dCoefficients.set(0, phiLeft.apply(t) + (u.get(k - 1).get(0) - 2. * u.get(k).get(0)) *
                        h * alpha / 2. / Math.pow(tau, 2));
                aCoefficients.set(N, -1. * gamma / h);
                bCoefficients.set(N, delta + gamma * (1. / h + 3. * h / 2. + h / 2. / Math.pow(tau, 2)));
                dCoefficients.set(N, phiRight.apply(t) + (2. * u.get(k).get(N) - u.get(k - 1).get(N)) *
                        h * gamma / 2. / Math.pow(tau, 2));
            } else {
                bCoefficients.set(0, beta - alpha / h);
                cCoefficients.set(0, alpha / h - alpha / 2. / h / sigma - 3. * h * alpha / 2.);
                dCoefficients.set(0, phiLeft.apply(t) + alpha / 2. / h / sigma *
                        (u.get(k - 1).get(1) - 2. * u.get(k).get(1)));
                aCoefficients.set(N, gamma / 2. / h / sigma + 3. * gamma * h / 2. - gamma / h);
                bCoefficients.set(N, gamma / h + delta);
                dCoefficients.set(N, phiRight.apply(t) - gamma / 2. / h / sigma *
                        (u.get(k - 1).get(N - 1) - 2 * u.get(k).get(N - 1)));
            }

            u.set(k + 1, tridiagonalAlgo(aCoefficients, bCoefficients, cCoefficients, dCoefficients));
        }

        Answer answer = new Answer(leftBound, rightBound, tau, h, T, u);
        return gson.toJson(answer);
    }

    class Answer {
        private final double leftBound;
        private final double rightBound;
        private final double tau;
        private final double h;
        private final double T;
        private final List<List<Double>> u;

        public Answer(double leftBound, double rightBound, double tau, double h, double t, List<List<Double>> u) {
            this.leftBound = leftBound;
            this.rightBound = rightBound;
            this.tau = tau;
            this.h = h;
            T = t;
            this.u = u;
        }
    }

    public enum ApproximationBoundaryConditionType {
        TWO_POINT_FIRST_DEGREE,
        TWO_POINT_SECOND_DEGREE,
        THREE_POINT_SECOND_DEGREE
    }

    public enum ApproximationInitialConditionType {
        FIRST_DEGREE,
        SECOND_DEGREE
    }
}
