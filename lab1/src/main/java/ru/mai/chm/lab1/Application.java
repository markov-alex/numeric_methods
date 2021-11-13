package ru.mai.chm.lab1;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.function.BiFunction;
import java.util.function.Function;

public class Application {
    public static void main(String[] args) throws IOException {
        Scanner sc = new Scanner(System.in);
        int N, K;
        double a, T;
        double b = 0, c = 0, alpha = 0, beta = 1, gamma = 0, delta = 1;
        double leftBound = 0, rightBound = Math.PI;
        double theta = 0.5;
        BiFunction<Double, Double, Double> f = (Double x, Double t) -> {
            return 0.;
        };
        ThreeArgFunction<Double, Double, Double, Double> analyticalSolution = (Double x, Double t, Double A) -> {
            return Math.exp(-1 * A * t) * Math.cos(x);
        };
        BiFunction<Double, Double, Double> boundaryCondition1 = (Double t, Double A) -> {
            return Math.exp(-1 * A * t);
        };
        BiFunction<Double, Double, Double> boundaryCondition2 = (Double t, Double A) -> {
            return -1 * Math.exp(-1 * A * t);
        };
        Function<Double, Double> initialCondition = (Double x) -> {
            return Math.cos(x);
        };
        Solver.ApproximationType approximationType = Solver.ApproximationType.TWO_POINT_FIRST_DEGREE;
        System.out.println("Введите N - число шагов по x");
        N = sc.nextInt();
        System.out.println("Введите K - число шагов по t");
        K = sc.nextInt();
        System.out.println("Введите коэффициент a");
        a = sc.nextDouble();
        System.out.println("Введите T - время");
        T = sc.nextDouble();
        System.out.println("Выберите тип аппроксимации:\n1. Двухточечная аппроксимация с первым порядком");
        System.out.println("2. Двухточечная аппроксимация со вторым порядком");
        System.out.println("3. Трехточечная аппроксимация со вторым порядком");
        int tmp = sc.nextInt();
        switch (tmp) {
            case 1:
                approximationType = Solver.ApproximationType.TWO_POINT_FIRST_DEGREE;
                break;
            case 2:
                approximationType = Solver.ApproximationType.TWO_POINT_SECOND_DEGREE;
                break;
            case 3:
                approximationType = Solver.ApproximationType.THREE_POINT_SECOND_DEGREE;
                break;
        }

        Solver solver = new Solver(leftBound, rightBound, a, N, K, T, b, c, alpha, beta, gamma, delta,
                f, analyticalSolution, boundaryCondition1, boundaryCondition2, initialCondition);
        File file = new File("results");
        file.mkdir();
        file = new File("results/explicitScheme.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.explicitScheme(approximationType));
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
        file = new File("results/implicitScheme.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.implicitScheme(approximationType));
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
        file = new File("results/crankNicolson.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.crankNicolson(approximationType, theta));
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
    }
}
