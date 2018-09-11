#pragma once

#include <QLabel>

#include "window.h"

class Window;

class View : public QLabel {

public:
	explicit View(Window *window, QWidget *parent = 0);
protected:
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
private:
	// Constants
	static constexpr double rotationSpeed = 180.0/500.0; // Degrees/pixel
	static constexpr double scrollSpeed = 1.01; // Factor per degrees turned (grows exponentially)

	Window *window;
	int prevX;
	int prevY;
};
