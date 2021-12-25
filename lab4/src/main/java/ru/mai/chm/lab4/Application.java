package ru.mai.chm.lab4;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.function.BiFunction;

public class Application {
    public static void main(String[] args) throws IOException {
        Scanner sc = new Scanner(System.in);
        int N, M, K;
        double a = 0., T = 0.;
        double lx = 0., rx = Math.PI / 4., ly = 0., ry = Math.log(2.);

        System.out.println("Введите N - число шагов по x");
        N = sc.nextInt();
        System.out.println("Введите M - число шагов по y");
        M = sc.nextInt();
        System.out.println("Введите K - число шагов по t");
        K = sc.nextInt();
        System.out.println("Введите коэффициент a");
        a = sc.nextDouble();
        System.out.println("Введите T - время");
        T = sc.nextDouble();

        double finalA = a;
        BiFunction<Double, Double, Double> phi1= (Double y, Double t) -> {
            return Math.cosh(y) * Math.exp(-3. * finalA * t);
        };
        BiFunction<Double, Double, Double> phi2 = (Double y, Double t) -> {
            return 0.;
        };
        double finalA1 = a;
        BiFunction<Double, Double, Double> phi3 = (Double x, Double t) -> {
            return Math.cos(2. * x) * Math.exp(-3. * finalA1 * t);
        };
        double finalA2 = a;
        BiFunction<Double, Double, Double> phi4 = (Double x, Double t) -> {
            return 5. / 4. * Math.cos(2. * x) * Math.exp(-3. * finalA2 * t);
        };
        BiFunction<Double, Double, Double> ksi = (Double x, Double y) -> {
            return Math.cos(2. * x) * Math.cosh(y);
        };
        double finalA3 = a;
        ThreeArgFunction<Double, Double, Double, Double> analyticalSolution = (Double x, Double y, Double t) -> {
            return Math.cos(2. * x) * Math.cosh(y) * Math.exp(-3. * finalA3 * t);
        };

        Solver solver = new Solver(lx, rx, ly, ry, N, M, K, T, a, phi1, phi2, phi3, phi4, ksi, analyticalSolution);
        File file = new File("results");
        file.mkdir();
        file = new File("results/alternatingDirections.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.alternatingDirections());
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
        file = new File("results/fractionalSteps.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.fractionalSteps());
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
        file = new File("results/analyticalSolution.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.answerAnalyticalSolution());
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
    }
}
