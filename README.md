# Задача Бензиностанции

Сорс кода в този проект е решение на задача "Бензиностанции" с имплементация на C#. Задачата е поставена на желаещите да се освободят от изпит в курса по "Основи на Програмирането", първи семестър, специалност ИКН в Унибит.

Подробна информация относно задачата, нейния анализ и теоретичните основи за нейното решение е поместена във файла GasStations.pdf.

## Как да тестваме кода

Имплементацията на решението е написана на C# използвайки .Net Core 3.1 технологията.

1. Отворете unibit-op-tsp.sln с Visual Studio
2. Направете build на целия сълюшън
3. Ако няма грешки изберете един от 2та варианта

### Вариант 1

* Стартирайте Program.cs в проекта GasStations
* На конзолата в терминала въведете тестови данни според спецификацията на задачата
* Проверете изходните данни за коректност

### Вариант 2

* Стартирайте Program.cs в проекта GasStationsDemo
* програмата автоматично зарежда тестовите данни
* изхода на конзолата дава допълнителна, интересна информация
* кода може да се редактира за да бъде зареден нов тест
* в проекта Utils, в клас TestData има допълнителни тестови данни и информация за верните им отховори.

4. Ако има грешки в build-а, моля отворете ново issue в сраницата на проекта в [GitHub](www.github.com/KostaDinkov/unibit-op-tsp).
