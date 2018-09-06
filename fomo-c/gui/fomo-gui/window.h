#pragma once

#include <functional>

#include <QMainWindow>
#include <QLabel>
#include <QFormLayout>
#include <QValidator>

#include "logger.h"
#include "FoMo.h"
#include "json.h"

using json = nlohmann::json;

namespace Ui {
	class MainWindow;
}

class Window : public QMainWindow {
	Q_OBJECT

public:
	explicit Window(QApplication &app, QWidget *parent = 0);
	~Window();

private:

	Ui::MainWindow *ui;
	QImage *image;
	FoMo::FoMoObject *FMO;

	// Managed by Qt, no need to call destructor
	Logger *log;
	QLabel *label;
	QTabWidget *controlPanel;
	QWidget *currentTab;
	QFormLayout *currentTabLayout;

	// Constants
	const std::string configPath = "../config.json";
	static constexpr double minAllowedValue = 1e-9; // Used for strictly positive values
	static constexpr double maxAllowedValue = 1e9;
	json config;
	double time = 0;

	// GUI creation methods
	void createTab(QString title);
	void addStringField(QString title, std::string name);
	void addIntegerField(QString title, std::string name, int min, int max = maxAllowedValue);
	void addDoubleField(QString title, std::string name, double min = minAllowedValue, double max = maxAllowedValue);
	void addButton(QString title, std::function<void (void)> handler);

private slots:
	void PrintTest();

};
