package ru.mai.chm.lab1;

@FunctionalInterface
public interface ThreeArgFunction<T1, T2, T3, R> {
    R apply(T1 t1, T2 t2, T3 t3);
}
