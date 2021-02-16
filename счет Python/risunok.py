# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.
"""
import matplotlib.pyplot as plt

#fig = plt.figure()
# Добавление на рисунок прямоугольной (по умолчанию) области рисования
scatter1 = plt.scatter(0.0, 1.0)
print('Scatter: ', type(scatter1))

graph1 = plt.plot([0.0, 0.0], [0.0, 1.0])
print('Plot: ', len(graph1), graph1)

text1 = plt.text(0.5, 0.5, 'Text on figure')
print('Text: ', type(text1))

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.show()