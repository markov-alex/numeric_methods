package ru.mai.chm.lab2;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Supplier;

public class Application {
    public static void main(String[] args) throws IOException {
        Scanner sc = new Scanner(System.in);
        int N, K;
        double T;
        double alpha = 0., beta = 1., gamma = 0., delta = 1.;
        double leftBound = 0., rightBound = Math.PI;
        BiFunction<Double, Double, Double> analyticalSolution = (Double x, Double t) -> {
            return Math.cos(x) * Math.sin(2. * t);
        };
        Function<Double, Double> phiLeft = (Double t) -> {
            return Math.sin(2. * t);
        };
        Function<Double, Double> phiRight = (Double t) -> {
            return -1. * Math.sin(2. * t);
        };
        Supplier<Double> psi1 = () -> {
            return 0.;
        };
        Function<Double, Double> psi2 = (Double x) -> {
            return 2. * Math.cos(x);
        };
        Supplier<Double> psi1SecondDerivative = () -> {
            return 0.;
        };
        Solver.ApproximationInitialConditionType approximationInitialConditionType =
                Solver.ApproximationInitialConditionType.FIRST_DEGREE;
        Solver.ApproximationBoundaryConditionType approximationBoundaryConditionType =
                Solver.ApproximationBoundaryConditionType.TWO_POINT_FIRST_DEGREE;
        System.out.println("Введите N - число шагов по x");
        N = sc.nextInt();
        System.out.println("Введите K - число шагов по t");
        K = sc.nextInt();
        System.out.println("Введите T - время");
        T = sc.nextDouble();
        System.out.println("Выберите тип аппроксимации первого слоя:");
        System.out.println("1. Первого порядка");
        System.out.println("2. Второго порядка");
        int tmp = sc.nextInt();
        switch (tmp) {
            case 1:
                approximationInitialConditionType = Solver.ApproximationInitialConditionType.FIRST_DEGREE;
                break;
            case 2:
                approximationInitialConditionType = Solver.ApproximationInitialConditionType.SECOND_DEGREE;
                break;
        }
        System.out.println("Выберите тип аппроксимации граничных условий:");
        System.out.println("1. Двухточечная аппроксимация с первым порядком");
        System.out.println("2. Двухточечная аппроксимация со вторым порядком");
        System.out.println("3. Трехточечная аппроксимация со вторым порядком");
        tmp = sc.nextInt();
        switch (tmp) {
            case 1:
                approximationBoundaryConditionType = Solver.ApproximationBoundaryConditionType.TWO_POINT_FIRST_DEGREE;
                break;
            case 2:
                approximationBoundaryConditionType = Solver.ApproximationBoundaryConditionType.TWO_POINT_SECOND_DEGREE;
                break;
            case 3:
                approximationBoundaryConditionType = Solver.ApproximationBoundaryConditionType.THREE_POINT_SECOND_DEGREE;
                break;
        }

        Solver solver = new Solver(leftBound, rightBound, N, K, T, alpha, beta, gamma, delta, analyticalSolution,
                phiLeft, phiRight, psi1, psi2, psi1SecondDerivative);
        File file = new File("results");
        file.mkdir();
        file = new File("results/explicitScheme.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.explicitScheme(approximationInitialConditionType,
                    approximationBoundaryConditionType));
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
        file = new File("results/implicitScheme.json");
        file.createNewFile();
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(solver.implicitScheme(approximationInitialConditionType,
                    approximationBoundaryConditionType));
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
