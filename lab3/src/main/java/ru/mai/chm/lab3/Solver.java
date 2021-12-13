package ru.mai.chm.lab3;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;

public class Solver {

    private static final Gson gson = new GsonBuilder()
            .setPrettyPrinting()
            .create();

    private final double lx;
    private final double rx;
    private final double ly;
    private final double ry;
    private final int N1;
    private final int N2;
    private final Function<Double, Double> f1;
    private final Function<Double, Double> f2;
    private final Function<Double, Double> f3;
    private final Function<Double, Double> f4;
    private final BiFunction<Double, Double, Double> analyticalSolution;

    public Solver(double lx, double rx, double ly, double ry, int N1, int N2,
                  Function<Double, Double> f1, Function<Double, Double> f2,
                  Function<Double, Double> f3, Function<Double, Double> f4,
                  BiFunction<Double, Double, Double> analyticalSolution) {
        this.lx = lx;
        this.rx = rx;
        this.ly = ly;
        this.ry = ry;
        this.N1 = N1;
        this.N2 = N2;
        this.f1 = f1;
        this.f2 = f2;
        this.f3 = f3;
        this.f4 = f4;
        this.analyticalSolution = analyticalSolution;
    }

    public String answerAnalyticalSolution() {
        double h1 = (rx - lx) / N1;
        double h2 = (ry - ly) / N2;

        List<List<Double>> u = new ArrayList<>(N1 + 1);
        for (int i = 0; i <= N1; ++i) {
            u.add(new ArrayList<>(Collections.nCopies(N2 + 1, 0.)));
        }

        for (int i = 0; i <= N1; ++i) {
            double x = i * h1 + lx;
            for (int j = 0; j <= N2; ++j) {
                double y = j * h2 + ly;
                u.get(i).set(j, analyticalSolution.apply(x, y));
            }
        }

        Answer answer = new Answer(lx, rx, ly, ry, h1, h2, 1, u);
        return gson.toJson(answer);
    }

    public List<List<Double>> initialization() {
        double h1 = (rx - lx) / N1;
        double h2 = (ry - ly) / N2;

        List<List<Double>> u = new ArrayList<>(N1 + 1);
        for (int i = 0; i <= N1; ++i) {
            u.add(new ArrayList<>(Collections.nCopies(N2 + 1, 0.)));
        }

        // левая граница
        for (int j = 0; j <= N2; ++j) {
            double y = j * h2 + ly;
            u.get(0).set(j, Math.cos(y));
        }
        // правая граница
        for (int j = 0; j <= N2; ++j) {
            double y = j * h2 + ly;
            u.get(N1).set(j, Math.E * Math.cos(y));
        }
        // внутренние точки
        for (int i = 1; i < N1; ++i) {
            double x = i * h1 + lx;
            for (int j = 1; j < N2; ++j) {
                double y = j * h2 + ly;
                u.get(i).set(j, (f2.apply(y) - f1.apply(y)) / (rx - lx) * (x - lx));
            }
        }
        // нижняя граница
        for (int i = 1; i < N1; ++i) {
            double x = i * h1 + lx;
            u.get(i).set(0, u.get(i).get(1) - h2 * f3.apply(x));
        }
        // верхняя граница
        for (int i = 1; i < N1; ++i) {
            double x = i * h1 + lx;
            u.get(i).set(N2, u.get(i).get(N2 - 1) + h2 * f4.apply(x));
        }

        return u;
    }

    private double norm(List<List<Double>> ukPlus1, List<List<Double>> uk) {
        double res = -1.;
        for (int i = 0; i <= N1; ++i) {
            for (int j = 0; j <= N2; ++j) {
                res = Math.max(Math.abs(ukPlus1.get(i).get(j) - uk.get(i).get(j)), res);
            }
        }
        return res;
    }

    private double norm2(List<Double> ukPlus1, List<Double> uk) {
        double res = -1.;
        for (int i = 0; i <= N1; ++i) {
            for (int j = 0; j <= N2; ++j) {
                res = Math.max(Math.abs(ukPlus1.get(mapCords(i, j)) - uk.get(mapCords(i, j))), res);
            }
        }
        return res;
    }

    private List<List<Double>> copy(List<List<Double>> u) {
        int n = u.size();
        List<List<Double>> res = new ArrayList<>(n);
        for (int i = 0; i < n; ++i) {
            List<Double> row = new ArrayList<>(u.get(i));
            res.add(row);
        }
        return res;
    }

    private List<Double> copy2(List<Double> u) {
        return new ArrayList<>(u);
    }

    public String liebmann(double eps) {
        double h1 = (rx - lx) / N1;
        double h2 = (ry - ly) / N2;
        List<List<Double>> uPrev = initialization();
        List<List<Double>> uCur = new ArrayList<>(N1 + 1);
        for (int i = 0; i <= N1; ++i) {
            uCur.add(new ArrayList<>(Collections.nCopies(N2 + 1, 0.)));
        }
        // левая граница
        for (int j = 0; j <= N2; ++j) {
            double y = j * h2 + ly;
            uCur.get(0).set(j, Math.cos(y));
        }
        // правая граница
        for (int j = 0; j <= N2; ++j) {
            double y = j * h2 + ly;
            uCur.get(N1).set(j, Math.E * Math.cos(y));
        }
        int k;
        for (k = 0; ; ++k) {
            // внутренние точки
            for (int i = 1; i < N1; ++i) {
                for (int j = 1; j < N2; ++j) {
                    uCur.get(i).set(j, (Math.pow(h2, 2) * uPrev.get(i + 1).get(j) +
                            Math.pow(h2, 2) * uPrev.get(i - 1).get(j) + Math.pow(h1, 2) * uPrev.get(i).get(j + 1) +
                            Math.pow(h1, 2) * uPrev.get(i).get(j - 1)) / 2 / (Math.pow(h1, 2) + Math.pow(h2, 2)));
                }
            }
            // нижняя линия
            for (int i = 1; i < N1; ++i) {
                double x = i * h1 + lx;
                uCur.get(i).set(0, uCur.get(i).get(1) - h2 * f3.apply(x));
            }
            // верхняя граница
            for (int i = 1; i < N1; ++i) {
                double x = i * h1 + lx;
                uCur.get(i).set(N2, uCur.get(i).get(N2 - 1) + h2 * f4.apply(x));
            }
            if (norm(uCur, uPrev) <= eps) {
                ++k;
                break;
            }

            uPrev = copy(uCur);
        }

        Answer answer = new Answer(lx, rx, ly, ry, h1, h2, k + 1, uCur);
        return gson.toJson(answer);
    }

    private int mapCords(int i, int j) {
        return i + j * (N1 + 1);
    }

    private int[] reverseMapCords(int c) {
        return new int[]{c % (N1 + 1), c / (N1 + 1)};
    }

    private List<Double> multMatrixVec(List<List<Double>> alpha, List<Double> u) {
        List<Double> res = new ArrayList<>(Collections.nCopies(u.size(), 0.));
        for (int i = 0; i < alpha.size(); ++i) {
            for (int j = 0; j < alpha.get(i).size(); ++j) {
                res.set(i, res.get(i) + alpha.get(i).get(j) * u.get(j));
            }
        }
        return res;
    }

    private List<Double> sumVec(List<Double> vec1, List<Double> vec2) {
        List<Double> res = new ArrayList<>(vec1);
        for (int i = 0; i < vec1.size(); ++i) {
            res.set(i, res.get(i) + vec2.get(i));
        }
        return res;
    }

    public String seidel(double eps) {
        double h1 = (rx - lx) / N1;
        double h2 = (ry - ly) / N2;
        List<List<Double>> uInit = initialization();

        List<Double> uPrev = new ArrayList<>(Collections.nCopies((N1 + 1) * (N2 + 1), 0.));
        List<Double> uCur = new ArrayList<>(Collections.nCopies((N1 + 1) * (N2 + 1), 0.));
        for (int i = 0; i <= N1; ++i) {
            for (int j = 0; j <= N2; ++j) {
                uPrev.set(mapCords(i, j), uInit.get(i).get(j));
            }
        }

        List<Double> beta = new ArrayList<>(Collections.nCopies((N1 + 1) * (N2 + 1), 0.));
        List<List<Double>> alpha = new ArrayList<>((N1 + 1) * (N2 + 1));
        for (int i = 0; i < (N1 + 1) * (N2 + 1); ++i) {
            alpha.add(new ArrayList<>(Collections.nCopies((N1 + 1) * (N2 + 1), 0.)));
        }

        // левая и правая границы
        for (int j = 0; j <= N2; ++j) {
            double y = j * h2 + ly;
            beta.set(mapCords(0, j), f1.apply(y));
            beta.set(mapCords(N1, j), f2.apply(y));
        }
        // нижняя и верхняя границы
        for (int i = 1; i < N1; ++i) {
            double x = i * h1 + lx;
            beta.set(mapCords(i, 0), -1. * h2 * f3.apply(x));
            beta.set(mapCords(i, N2), h2 * f4.apply(x));
            alpha.get(mapCords(i, 0)).set(mapCords(i, 1), 1.);
            alpha.get(mapCords(i, N2)).set(mapCords(i, N2 - 1), 1.);
        }
        // внутренние точки
        for (int i = 1; i < N1; ++i) {
            for (int j = 1; j < N2; ++j) {
                alpha.get(mapCords(i, j)).
                        set(mapCords(i + 1, j), Math.pow(h2, 2) / 2. / (Math.pow(h1, 2) + Math.pow(h2, 2)));
                alpha.get(mapCords(i, j)).
                        set(mapCords(i - 1, j), Math.pow(h2, 2) / 2. / (Math.pow(h1, 2) + Math.pow(h2, 2)));
                alpha.get(mapCords(i, j)).
                        set(mapCords(i, j + 1), Math.pow(h1, 2) / 2. / (Math.pow(h1, 2) + Math.pow(h2, 2)));
                alpha.get(mapCords(i, j)).
                        set(mapCords(i, j - 1), Math.pow(h1, 2) / 2. / (Math.pow(h1, 2) + Math.pow(h2, 2)));
            }
        }

        int k;
        for (k = 0; ; ++k) {
            List<Double> tmpU = copy2(uPrev);
            for (int q = 0; q < tmpU.size(); ++q) {
                List<Double> newU = multMatrixVec(alpha, tmpU);
                tmpU.set(q, sumVec(newU, beta).get(q));
            }
            uCur = tmpU;

            if (norm2(uCur, uPrev) <= eps) {
                ++k;
                break;
            }

            uPrev = copy2(uCur);
        }

        List<List<Double>> u = new ArrayList<>(N1 + 1);
        for (int i = 0; i <= N1; ++i) {
            u.add(new ArrayList<>(Collections.nCopies(N2 + 1, 0.)));
        }
        for (int q = 0; q < uCur.size(); ++q) {
            int i = reverseMapCords(q)[0], j = reverseMapCords(q)[1];
            u.get(i).set(j, uCur.get(q));
        }

        Answer answer = new Answer(lx, rx, ly, ry, h1, h2, k + 1, u);
        return gson.toJson(answer);
    }

    public String relaxation(double eps, double tau) {
        double h1 = (rx - lx) / N1;
        double h2 = (ry - ly) / N2;
        List<List<Double>> uInit = initialization();

        List<Double> uPrev = new ArrayList<>(Collections.nCopies((N1 + 1) * (N2 + 1), 0.));
        List<Double> uCur = new ArrayList<>(Collections.nCopies((N1 + 1) * (N2 + 1), 0.));
        for (int i = 0; i <= N1; ++i) {
            for (int j = 0; j <= N2; ++j) {
                uPrev.set(mapCords(i, j), uInit.get(i).get(j));
            }
        }

        List<Double> beta = new ArrayList<>(Collections.nCopies((N1 + 1) * (N2 + 1), 0.));
        List<List<Double>> alpha = new ArrayList<>((N1 + 1) * (N2 + 1));
        for (int i = 0; i < (N1 + 1) * (N2 + 1); ++i) {
            alpha.add(new ArrayList<>(Collections.nCopies((N1 + 1) * (N2 + 1), 0.)));
        }

        // левая и правая границы
        for (int j = 0; j <= N2; ++j) {
            double y = j * h2 + ly;
            beta.set(mapCords(0, j), f1.apply(y));
            beta.set(mapCords(N1, j), f2.apply(y));
        }
        // нижняя и верхняя границы
        for (int i = 1; i < N1; ++i) {
            double x = i * h1 + lx;
            beta.set(mapCords(i, 0), -1. * h2 * f3.apply(x));
            beta.set(mapCords(i, N2), h2 * f4.apply(x));
            alpha.get(mapCords(i, 0)).set(mapCords(i, 1), 1.);
            alpha.get(mapCords(i, N2)).set(mapCords(i, N2 - 1), 1.);
        }
        // внутренние точки
        for (int i = 1; i < N1; ++i) {
            for (int j = 1; j < N2; ++j) {
                alpha.get(mapCords(i, j)).
                        set(mapCords(i + 1, j), tau * Math.pow(h2, 2) / 2. / (Math.pow(h1, 2) + Math.pow(h2, 2)));
                alpha.get(mapCords(i, j)).
                        set(mapCords(i - 1, j), tau * Math.pow(h2, 2) / 2. / (Math.pow(h1, 2) + Math.pow(h2, 2)));
                alpha.get(mapCords(i, j)).
                        set(mapCords(i, j + 1), tau * Math.pow(h1, 2) / 2. / (Math.pow(h1, 2) + Math.pow(h2, 2)));
                alpha.get(mapCords(i, j)).
                        set(mapCords(i, j - 1), tau * Math.pow(h1, 2) / 2. / (Math.pow(h1, 2) + Math.pow(h2, 2)));
            }
        }

        int k;
        for (k = 0; ; ++k) {
            List<Double> tmpU = copy2(uPrev);
            for (int q = 0; q < tmpU.size(); ++q) {
                int[] cords = reverseMapCords(q);
                if (cords[0] >= 1 && cords[0] <= N1 - 1 && cords[1] >= 1 && cords[1] <= N2 - 1) {
                    beta.set(q, (1 - tau) * tmpU.get(q));
                }
            }
            for (int q = 0; q < tmpU.size(); ++q) {
                List<Double> newU = multMatrixVec(alpha, tmpU);
                tmpU.set(q, sumVec(newU, beta).get(q));
            }
            uCur = tmpU;

            if (norm2(uCur, uPrev) <= eps) {
                ++k;
                break;
            }

            uPrev = copy2(uCur);
        }

        List<List<Double>> u = new ArrayList<>(N1 + 1);
        for (int i = 0; i <= N1; ++i) {
            u.add(new ArrayList<>(Collections.nCopies(N2 + 1, 0.)));
        }
        for (int q = 0; q < uCur.size(); ++q) {
            int i = reverseMapCords(q)[0], j = reverseMapCords(q)[1];
            u.get(i).set(j, uCur.get(q));
        }

        Answer answer = new Answer(lx, rx, ly, ry, h1, h2, k + 1, u);
        return gson.toJson(answer);
    }

    class Answer {
        private final double lx;
        private final double rx;
        private final double ly;
        private final double ry;
        private final double h1;
        private final double h2;
        private final int cntIter;
        private final List<List<Double>> u;

        public Answer(double lx, double rx, double ly, double ry, double h1, double h2, int cntIter,
                      List<List<Double>> u) {
            this.lx = lx;
            this.rx = rx;
            this.ly = ly;
            this.ry = ry;
            this.h1 = h1;
            this.h2 = h2;
            this.cntIter = cntIter;
            this.u = u;
        }
    }
}
