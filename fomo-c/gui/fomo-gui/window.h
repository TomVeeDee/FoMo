#pragma once

#include <functional>
#include <chrono>

#include <QMainWindow>
#include <QPainter>
#include <QLabel>
#include <QFormLayout>
#include <QLineEdit>

#include "logger.h"
#include "view.h"
#include "FoMo.h"
#include "json.h"

using json = nlohmann::json;

class View;

class Window : public QMainWindow {
	Q_OBJECT

public:

	explicit Window(QApplication &app, QWidget *parent = 0);
	~Window();

	// View parameter updating methods
	void updateAngle(double dl, double db);
	void updateViewSize(double factor);

private:

	// Constants
	const std::string configPath = "config.json";
	static constexpr double minAllowedValue = 1e-9; // Used for strictly positive values
	static constexpr double maxAllowedValue = 1e9;
	static constexpr int views_amount = 3;

	// General
	std::chrono::time_point<std::chrono::high_resolution_clock> start;
	json config;
	bool updatedRenderingSettings = true;
	bool updatedViewParameters = true;
	double renderTimeLastSecond = 0;
	int framesLastSecond = 0;

	// State variables
	int rendererState = 0;

	// Need to be destroyed manually
	FoMo::FoMoObject *FMO = NULL; // Stores the GoftCube without making further copies
	FoMo::RegularGridRendererWrapper *renderer = NULL;
	unsigned char *imageBuffer = NULL;

	// Managed by Qt, no need to call destructor
	Logger *logger;
	View *views[views_amount];
	QTimer *renderTimer;
	QTabWidget *controlPanel;
	QWidget *currentTab;
	QFormLayout *currentTabLayout;
	QLineEdit *l_field;
	QLineEdit *b_field;
	QLineEdit *view_width_field;
	QLineEdit *view_height_field;

	// GUI creation methods
	void createTab(QString title);
	QLineEdit* addStringField(QString title, std::string name);
	QLineEdit* addIntegerField(QString title, std::string name, int min, int max = maxAllowedValue);
	QLineEdit* addDoubleField(QString title, std::string name, double min = minAllowedValue, double max = maxAllowedValue);
	void addButton(QString title, std::function<void (void)> handler);

	// Handlers
	void readAndPreprocessDataCube();
	void constructRegularGrid();
	void updateRenderingSettings();
	void renderToFile();
	void renderToView();

	// Helper methods
	inline std::chrono::time_point<std::chrono::high_resolution_clock> time_now();
	inline float parseFloat(std::ifstream &in, char *buffer, const int bufferSize, int &bufferPos, int &charsRead);
	inline float toRadians(float degrees);
	inline void log(const std::string &text);
	inline void logProcessStart();
	inline void logProcessFinished(std::string name);

};
