# Splines


В проекте теперь реализован градиентный спуск с методом наискорейшего спуска(заключается в поиске оптимального шага, чтобы была сходимость и выполнялось минимальное число итераций)
В целом, можно было обойтись и без этого метода, просто взяв очень маленькое eps, но тогда количество итераций было бы намного больше.

Программе нужно ввести два файла, в них лежат пары координат точек для каждого сплайна. Программа также создает файл с точками, чтобы впоследствии их можно было нарисовать
в питоне(Питоновский файл надо положить в одну папку с exe чтобы работало). Для примера в файле лежат файлы 1.txt и 2.txt. В них пример непересекающихся сплайнов, 
в которых ищется расстояние
