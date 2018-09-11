#include <iostream>
#include <math.h>

#include <QMouseEvent>
#include <QWheelEvent>
#include <QPoint>

#include "view.h"
#include "window.h"

View::View(Window *window, QWidget *parent) : QLabel(parent) {
	this->window = window;
}

void View::mousePressEvent(QMouseEvent *event) {
	prevX = event->pos().x();
	prevY = event->pos().y();
}

void View::mouseMoveEvent(QMouseEvent *event) {
	int x = event->pos().x();
	int y = event->pos().y();
	window->updateAngle((y - prevY)*rotationSpeed, (x - prevX)*rotationSpeed);
	prevX = x;
	prevY = y;
}

void View::wheelEvent(QWheelEvent *event) {
	window->updateViewSize(pow(scrollSpeed, -event->delta()/8.0));
}
