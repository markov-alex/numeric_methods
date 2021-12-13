package ru.mai.chm.lab3;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.function.BiFunction;
import java.util.function.Function;

public class Application {
    public static void main(String[] args) throws IOException {
        Scanner sc = new Scanner(System.in);
        double lx = 0.;
        double rx = 1.;
        double ly = 0.;
        double ry = Math.PI / 2.;
        int N1;
        int N2;
        Function<Double, Double> f1 = (Double y) -> {
            return Math.cos(y);
        };
        Function<Double, Double> f2 = (Double y) -> {
            return Math.E * Math.cos(y);
        };
        Function<Double, Double> f3 = (Double x) -> {
            return 0.;
        };
        Function<Double, Double> f4 = (Double x) -> {
            return -1. * Math.exp(x);
        };
        BiFunction<Double, Double, Double> analyticalSolution = (Double x, Double y) -> {
            return Math.exp(x) * Math.cos(y);
        };
        double eps = 0.001;
        double tau = 0.8;
        System.out.println("Введите N1 - число шагов по x");
        N1 = sc.nextInt();
        System.out.println("Введите N2 - число шагов по y");
        N2 = sc.nextInt();
        System.out.println("Введите epsilon");
        eps = sc.nextDouble();
        System.out.println("Введите tau");
        tau = sc.nextDouble();
        Solver solver = new Solver(lx, rx, ly, ry, N1, N2, f1, f2, f3, f4, analyticalSolution);

        File file = new File("results");
        file.mkdir();
        file = new File("results/analyticalSolution.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.answerAnalyticalSolution());
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
        file = new File("results/liebmann.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.liebmann(eps));
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
        file = new File("results/seidel.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.seidel(eps));
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
        file = new File("results/relaxation.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.relaxation(eps, tau));
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
    }
}
