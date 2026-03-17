# Модель движения колеса по деформируемому основанию на C++
Модель синтезирована на основе [статьи](https://cyberleninka.ru/article/n/matematicheskaya-model-silovogo-vzaimodeystviya-kolesa-s-gruntom-pri-povorote-mashiny), посвящённой силовому взаимодействию колеса с грунтом в пятне контакта, и [статьи](https://cyberleninka.ru/article/n/matematicheskaya-model-pryamolineynogo-kacheniya-elastichnogo-kolesa-po-nerovnomu-deformiruemomu-opornomu-osnovaniyu) с моделированием качения по деформируемому основанию.
Расчётная схема модели колеса, по которой писалась реализация, представлена на рисунке:

<img width="732" height="660" alt="image" src="https://github.com/user-attachments/assets/371a2662-5dd3-4a15-a5f2-7057c32400fe" />

> Для каждого момента времени формируется и решается система ОДУ (метод Рунге-Кутты 4-го порядка) на основе интеграла по пятну контакта (составная формула трапеций).
