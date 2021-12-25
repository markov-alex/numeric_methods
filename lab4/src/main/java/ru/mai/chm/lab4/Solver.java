package ru.mai.chm.lab4;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;

public class Solver {

    private static final Gson gson = new GsonBuilder()
            .setPrettyPrinting()
            .create();

    private final double lx;
    private final double rx;
    private final double ly;
    private final double ry;
    private final int N;
    private final int M;
    private final int K;
    private final double T;
    private final double a;
    private final BiFunction<Double, Double, Double> phi1;
    private final BiFunction<Double, Double, Double> phi2;
    private final BiFunction<Double, Double, Double> phi3;
    private final BiFunction<Double, Double, Double> phi4;
    private final BiFunction<Double, Double, Double> ksi;
    private final ThreeArgFunction<Double, Double, Double, Double> analyticalSolution;

    public Solver(double lx, double rx, double ly, double ry, int N, int M, int K, double T, double a,
                  BiFunction<Double, Double, Double> phi1, BiFunction<Double, Double, Double> phi2,
                  BiFunction<Double, Double, Double> phi3, BiFunction<Double, Double, Double> phi4,
                  BiFunction<Double, Double, Double> ksi,
                  ThreeArgFunction<Double, Double, Double, Double> analyticalSolution) {
        this.lx = lx;
        this.rx = rx;
        this.ly = ly;
        this.ry = ry;
        this.N = N;
        this.M = M;
        this.K = K;
        this.T = T;
        this.a = a;
        this.phi1 = phi1;
        this.phi2 = phi2;
        this.phi3 = phi3;
        this.phi4 = phi4;
        this.ksi = ksi;
        this.analyticalSolution = analyticalSolution;
    }

    public String answerAnalyticalSolution() {
        double tau = T / K;
        double h1 = (rx - lx) / N;
        double h2 = (ry - ly) / M;

        List<List<List<Double>>> u = new ArrayList<>(K + 1);
        for (int k = 0; k <= K; ++k) {
            List<List<Double>> tmp = new ArrayList<>(N + 1);
            for (int i = 0; i <= N; ++i) {
                tmp.add(new ArrayList<>(Collections.nCopies(M + 1, 0.)));
            }
            u.add(tmp);
        }

        double x, y, t;
        for (int k = 0; k <= K; ++k) {
            t = k * tau;
            for (int i = 0; i <= N; ++i) {
                x = lx + i * h1;
                for (int j = 0; j <= M; ++j) {
                    y = ly + j * h2;
                    u.get(k).get(i).set(j, analyticalSolution.apply(x, y, t));
                }
            }
        }

        Answer answer = new Answer(lx, rx, ly, ry, T, tau, h1, h2, u);
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

    public String alternatingDirections() {
        double tau = T / K;
        double h1 = (rx - lx) / N;
        double h2 = (ry - ly) / M;

        List<List<List<Double>>> u = new ArrayList<>(K + 1);
        for (int k = 0; k <= K; ++k) {
            List<List<Double>> tmp = new ArrayList<>(N + 1);
            for (int i = 0; i <= N; ++i) {
                tmp.add(new ArrayList<>(Collections.nCopies(M + 1, 0.)));
            }
            u.add(tmp);
        }

        // t = 0
        for (int i = 0; i <= N; ++i) {
            double x = lx + i * h1;
            for (int j = 0; j <= M; ++j) {
                double y = ly + j * h2;
                u.get(0).get(i).set(j, ksi.apply(x, y));
            }
        }

        for (int k = 0; k < K; ++k) {
            double t = (k + 1) * tau;
            List<List<Double>> kHalf = new ArrayList<>(N + 1);
            for (int i = 0; i <= N; ++i) {
                kHalf.add(new ArrayList<>(Collections.nCopies(M + 1, 0.)));
            }

            for (int j = 1; j <= M - 1; ++j) {
                double y = ly + j * h2;
                List<Double> aCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));
                List<Double> bCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));
                List<Double> cCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));
                List<Double> dCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));

                for (int i = 1; i <= N - 1; ++i) {
                    aCoefficients.set(i, a / Math.pow(h1, 2));
                    bCoefficients.set(i, -2. * a / Math.pow(h1, 2) - 2. / tau);
                    cCoefficients.set(i, a / Math.pow(h1, 2));
                    dCoefficients.set(i, a / Math.pow(h2, 2) * (-1. * u.get(k).get(i).get(j + 1)
                            + 2 * u.get(k).get(i).get(j) - u.get(k).get(i).get(j - 1))
                            - 2. / tau * u.get(k).get(i).get(j));
                }

                aCoefficients.set(0, 0.);
                bCoefficients.set(0, 1.);
                cCoefficients.set(0, 0.);
                dCoefficients.set(0, phi1.apply(y, t - tau / 2.));
                aCoefficients.set(N, 0.);
                bCoefficients.set(N, 1.);
                cCoefficients.set(N, 0.);
                dCoefficients.set(N, phi2.apply(y, t - tau / 2.));

                List<Double> ans = tridiagonalAlgo(aCoefficients, bCoefficients, cCoefficients, dCoefficients);
                for (int i = 0; i <= N; ++i) {
                    kHalf.get(i).set(j, ans.get(i));
                }
            }

            for (int i = 0; i <= N; ++i) {
                double x = lx + i * h1;
                kHalf.get(i).set(0, phi3.apply(x, t - tau / 2.));
                kHalf.get(i).set(M, phi4.apply(x, t - tau / 2.));
            }

            for (int i = 1; i <= N - 1; ++i) {
                double x = lx + i * h1;
                List<Double> aCoefficients = new ArrayList<>(Collections.nCopies(M + 1, 0.));
                List<Double> bCoefficients = new ArrayList<>(Collections.nCopies(M + 1, 0.));
                List<Double> cCoefficients = new ArrayList<>(Collections.nCopies(M + 1, 0.));
                List<Double> dCoefficients = new ArrayList<>(Collections.nCopies(M + 1, 0.));

                for (int j = 1; j <= M - 1; ++j) {
                    aCoefficients.set(j, a / Math.pow(h2, 2));
                    bCoefficients.set(j, -2. * a / Math.pow(h2, 2) - 2. / tau);
                    cCoefficients.set(j, a / Math.pow(h2, 2));
                    dCoefficients.set(j, a / Math.pow(h1, 2)
                            * (-1. * kHalf.get(i + 1).get(j) + 2 * kHalf.get(i).get(j) - kHalf.get(i - 1).get(j))
                            - 2. / tau * kHalf.get(i).get(j));
                }

                aCoefficients.set(0, 0.);
                bCoefficients.set(0, 1.);
                cCoefficients.set(0, 0.);
                dCoefficients.set(0, phi3.apply(x, t));
                aCoefficients.set(M, 0.);
                bCoefficients.set(M, 1.);
                cCoefficients.set(M, 0.);
                dCoefficients.set(M, phi4.apply(x, t));

                List<Double> ans = tridiagonalAlgo(aCoefficients, bCoefficients, cCoefficients, dCoefficients);
                for (int j = 0; j <= M; ++j) {
                    u.get(k + 1).get(i).set(j, ans.get(j));
                }
            }

            for (int j = 0; j <= M; ++j) {
                double y = ly + j * h2;
                u.get(k + 1).get(0).set(j, phi1.apply(y, t));
                u.get(k + 1).get(N).set(j, phi2.apply(y, t));
            }
        }

        Answer answer = new Answer(lx, rx, ly, ry, T, tau, h1, h2, u);
        return gson.toJson(answer);
    }

    public String fractionalSteps() {
        double tau = T / K;
        double h1 = (rx - lx) / N;
        double h2 = (ry - ly) / M;

        List<List<List<Double>>> u = new ArrayList<>(K + 1);
        for (int k = 0; k <= K; ++k) {
            List<List<Double>> tmp = new ArrayList<>(N + 1);
            for (int i = 0; i <= N; ++i) {
                tmp.add(new ArrayList<>(Collections.nCopies(M + 1, 0.)));
            }
            u.add(tmp);
        }

        // t = 0
        for (int i = 0; i <= N; ++i) {
            double x = lx + i * h1;
            for (int j = 0; j <= M; ++j) {
                double y = ly + j * h2;
                u.get(0).get(i).set(j, ksi.apply(x, y));
            }
        }

        for (int k = 0; k < K; ++k) {
            double t = (k + 1) * tau;
            List<List<Double>> kHalf = new ArrayList<>(N + 1);
            for (int i = 0; i <= N; ++i) {
                kHalf.add(new ArrayList<>(Collections.nCopies(M + 1, 0.)));
            }

            for (int j = 1; j <= M - 1; ++j) {
                double y = ly + j * h2;
                List<Double> aCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));
                List<Double> bCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));
                List<Double> cCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));
                List<Double> dCoefficients = new ArrayList<>(Collections.nCopies(N + 1, 0.));

                for (int i = 1; i <= N - 1; ++i) {
                    aCoefficients.set(i, a / Math.pow(h1, 2));
                    bCoefficients.set(i, -2. * a / Math.pow(h1, 2) - 1 / tau);
                    cCoefficients.set(i, a / Math.pow(h1, 2));
                    dCoefficients.set(i, u.get(k).get(i).get(j) * -1. / tau);
                }

                aCoefficients.set(0, 0.);
                bCoefficients.set(0, 1.);
                cCoefficients.set(0, 0.);
                dCoefficients.set(0, phi1.apply(y, t - tau / 2.));
                aCoefficients.set(N, 0.);
                bCoefficients.set(N, 1.);
                cCoefficients.set(N, 0.);
                dCoefficients.set(N, phi2.apply(y, t - tau / 2.));

                List<Double> ans = tridiagonalAlgo(aCoefficients, bCoefficients, cCoefficients, dCoefficients);
                for (int i = 0; i <= N; ++i) {
                    kHalf.get(i).set(j, ans.get(i));
                }
            }

            for (int i = 0; i <= N; ++i) {
                double x = lx + i * h1;
                kHalf.get(i).set(0, phi3.apply(x, t - tau / 2.));
                kHalf.get(i).set(M, phi4.apply(x, t - tau / 2.));
            }

            for (int i = 1; i <= N - 1; ++i) {
                double x = lx + i * h1;
                List<Double> aCoefficients = new ArrayList<>(Collections.nCopies(M + 1, 0.));
                List<Double> bCoefficients = new ArrayList<>(Collections.nCopies(M + 1, 0.));
                List<Double> cCoefficients = new ArrayList<>(Collections.nCopies(M + 1, 0.));
                List<Double> dCoefficients = new ArrayList<>(Collections.nCopies(M + 1, 0.));

                for (int j = 1; j <= M - 1; ++j) {
                    aCoefficients.set(j, a / Math.pow(h2, 2));
                    bCoefficients.set(j, -2. * a / Math.pow(h2, 2) - 1. / tau);
                    cCoefficients.set(j, a / Math.pow(h2, 2));
                    dCoefficients.set(j, kHalf.get(i).get(j) * -1. / tau);
                }

                aCoefficients.set(0, 0.);
                bCoefficients.set(0, 1.);
                cCoefficients.set(0, 0.);
                dCoefficients.set(0, phi3.apply(x, t));
                aCoefficients.set(M, 0.);
                bCoefficients.set(M, 1.);
                cCoefficients.set(M, 0.);
                dCoefficients.set(M, phi4.apply(x, t));

                List<Double> ans = tridiagonalAlgo(aCoefficients, bCoefficients, cCoefficients, dCoefficients);
                for (int j = 0; j <= M; ++j) {
                    u.get(k + 1).get(i).set(j, ans.get(j));
                }
            }

            for (int j = 0; j <= M; ++j) {
                double y = ly + j * h2;
                u.get(k + 1).get(0).set(j, phi1.apply(y, t));
                u.get(k + 1).get(N).set(j, phi2.apply(y, t));
            }
        }

        Answer answer = new Answer(lx, rx, ly, ry, T, tau, h1, h2, u);
        return gson.toJson(answer);
    }

    class Answer {
        private final double lx;
        private final double rx;
        private final double ly;
        private final double ry;
        private final double T;
        private final double tau;
        private final double h1;
        private final double h2;
        private final List<List<List<Double>>> u;

        public Answer(double lx, double rx, double ly, double ry, double T, double tau, double h1, double h2,
                      List<List<List<Double>>> u) {
            this.lx = lx;
            this.rx = rx;
            this.ly = ly;
            this.ry = ry;
            this.T = T;
            this.tau = tau;
            this.h1 = h1;
            this.h2 = h2;
            this.u = u;
        }
    }
}
