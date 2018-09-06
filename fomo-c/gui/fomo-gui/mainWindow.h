#pragma once

#include <QMainWindow>
#include <QLabel>

#include "logger.h"
using json = nlohmann::json;
#include "FoMo.h"
#include "json.h"

namespace Ui {
	class MainWindow;
}

class MainWindow : public QMainWindow {
	Q_OBJECT

public:
	explicit MainWindow(QApplication &app, QWidget *parent = 0);
	~MainWindow();

public slots:
	void PrintTest();

private:

	Ui::MainWindow *ui;
	QImage *image;
	FoMo::FoMoObject *FMO;

	// Managed by Qt, no need to call destructor
	Logger *log;
	QLabel *label;

	static const std::string configPath = "config.json";
	json config;
	double time = 0;

};
