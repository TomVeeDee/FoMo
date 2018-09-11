#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <functional>
#include <cmath>
#include <algorithm>

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>

#include <QApplication>
#include <QStyle>
#include <QDesktopWidget>
#include <QGridLayout>
#include <QFormLayout>
#include <QTimer>
#include <QImage>
#include <QLabel>
#include <QLineEdit>
#include <QIntValidator>
#include <QDoubleValidator>
#include <QPushButton>

#include "window.h"
#include "ui_mainwindow.h"
#include "FoMo.h"

#define WINDOW_RENDERER_STATE_HAS_DATACUBE 1
#define WINDOW_RENDERER_STATE_HAS_REGULAR_GRID 2
#define WINDOW_RENDERER_STATE_HAS_RENDERING_SETTINGS 3

Window::Window(QApplication &app, QWidget *parent) : QMainWindow(parent) {

	// Initializes this window
	// Caller is still responsible for calling show() when necessary

	// General
	resize(1280, 720);
	setWindowTitle("FoMo rendering GUI");
	setGeometry(QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(), app.desktop()->availableGeometry())); // Centers window on screen

	// Load config
	// Values are sanity checked by the same methods that add them to the lay-out
	std::ifstream in(configPath);
	in >> config;
	in.close();

	// Top-level layout
	QWidget *centralWidget = new QWidget(this);
	QGridLayout *centralLayout = new QGridLayout(); // No need to assign parent here, setLayout already assigns ownership
	centralWidget->setLayout(centralLayout);
	setCentralWidget(centralWidget);

	// Control panel initialization
	controlPanel = new QTabWidget(centralWidget);
	centralLayout->addWidget(controlPanel, 0, 1);

	// Datacube
	createTab("DataCube");
	addStringField("Datacube file", "datacube_file");
	addStringField("CHIANTI file", "chianti_file");
	addStringField("Abundance file", "abundance_file");
	addIntegerField("Dimensions", "dims", 2, 3);
	addIntegerField("Variables", "nvars", 1);
	addIntegerField("X points", "nx", 2);
	addIntegerField("Y points", "ny", 2);
	addIntegerField("Z points", "nz", 1);
	addIntegerField("Data points", "ng", 1);
	addButton("Load and preprocess DataCube", [this] () {
		readAndPreprocessDataCube();
	});

	// Regular grid
	createTab("Regular grid");
	addIntegerField("X points", "gridx", 2);
	addIntegerField("Y points", "gridy", 2);
	addIntegerField("Z points", "gridz", 1);
	addDoubleField("Max X distance", "max_distance_x");
	addDoubleField("Max Y distance", "max_distance_y");
	addDoubleField("Max Z distance", "max_distance_z");
	addButton("Construct regular grid", [this] () {
		constructRegularGrid();
	});

	// Rendering settings
	createTab("Rendering settings");
	addIntegerField("X pixels", "x_pixel", 2);
	addIntegerField("Y pixels", "y_pixel", 2);
	addIntegerField("Lambda pixels", "lambda_pixel", 1);
	addDoubleField("Lambda width (m/s)", "lambda_width");
	addDoubleField("Maximal intensity (ergs/(cm^2*s*sr))", "max_intensity");
	addDoubleField("Target FPS", "fps");
	addButton("Update rendering settings", [this] () {
		updateRenderingSettings();
	});

	// Render
	createTab("Render");
	l_field = addDoubleField("l (degrees)", "l", -maxAllowedValue);
	b_field = addDoubleField("b (degrees)", "b", -maxAllowedValue);
	view_width_field = addDoubleField("View width (Mm)", "view_width");
	view_height_field = addDoubleField("View height (Mm)", "view_height");
	addButton("Update view parameters", [this] () {
		// Will re-render the image next time the rendering timer triggers
		updatedViewParameters = true;
	});
	addStringField("Output file", "outfile");
	addButton("Render to file", [this] () {
		renderToFile();
	});

	// Save config button
	QPushButton *save_config_button = new QPushButton("Save configuration to file", centralWidget);
	connect(save_config_button, &QPushButton::clicked, [this] () {
		std::ofstream out(this->configPath);
		out << std::setw(4) << this->config;
		out.close();
	});
	centralLayout->addWidget(save_config_button, 1, 1);

	// Log
	logger = new Logger(centralWidget);
	centralLayout->addWidget(logger, 2, 1);

	// View
	image = QImage(512, 512, QImage::Format_ARGB32);
	view = new View(this, centralWidget);
	centralLayout->addWidget(view, 0, 0, -1, 1, Qt::Alignment(Qt::AlignCenter));

	// Set up render timer
	QTimer *renderTimer = new QTimer(view);
	renderTimer->start(16);
	connect(renderTimer, &QTimer::timeout, this, [this] () {
		renderToView();
	});

}

Window::~Window() {
	// Most destruction is handled by the Qt parent-child hierarchy
	if (FMO != NULL)
		delete FMO;
	if (renderer != NULL)
		delete renderer;
	if (imageBuffer != NULL)
		delete[] imageBuffer;
}

// GUI creation methods

void Window::createTab(QString title) {
	// Creates a new tab with the given title
	currentTab = new QWidget(controlPanel);
	controlPanel->addTab(currentTab, title);
	currentTabLayout = new QFormLayout(); // No need to assign parent here, setLayout already assigns ownership
	currentTab->setLayout(currentTabLayout);
}

QLineEdit* Window::addStringField(QString title, std::string name) {
	// Adds a labelled string field connected to the config to the current tab
	QLineEdit *lineEdit = new QLineEdit(QString::fromStdString(static_cast<std::string const&>(config[name])), currentTab);
	connect(lineEdit, &QLineEdit::textEdited, [this, name] (const QString &text) {
		this->config[name] = text.toStdString();
	});
	currentTabLayout->addRow(title, lineEdit);
	return lineEdit;
}

QLineEdit* Window::addIntegerField(QString title, std::string name, int min, int max) {
	// Adds a labelled integer field connected to the config to the current tab
	// max is inclusive
	// Also does sanity check
	config[name] = std::max(min, std::min(max, static_cast<int>(config[name])));
	QLineEdit *lineEdit = new QLineEdit(QString::number(static_cast<int>(config[name])), currentTab);
	lineEdit->setValidator(new QIntValidator(min, max, lineEdit));
	connect(lineEdit, &QLineEdit::textEdited, [this, name, lineEdit] (const QString &text) {
		if (lineEdit->hasAcceptableInput())
			this->config[name] = text.toInt();
	});
	currentTabLayout->addRow(title, lineEdit);
	return lineEdit;
}

QLineEdit* Window::addDoubleField(QString title, std::string name, double min, double max) {
	// Adds a labelled double field connected to the config to the current tab
	// max is inclusive
	// Also does sanity check
	config[name] = std::max(min, std::min(max, static_cast<double>(config[name])));
	QLineEdit *lineEdit = new QLineEdit(QString::number(static_cast<double>(config[name])), currentTab);
	lineEdit->setValidator(new QDoubleValidator(min, max, 10, lineEdit));
	connect(lineEdit, &QLineEdit::textEdited, [this, name, lineEdit] (const QString &text) {
		if (lineEdit->hasAcceptableInput())
			this->config[name] = text.toDouble();
	});
	currentTabLayout->addRow(title, lineEdit);
	return lineEdit;
}

void Window::addButton(QString title, std::function<void (void)> handler) {
	// Adds a button to the current tab
	QPushButton *button = new QPushButton(title, currentTab);
	connect(button, &QPushButton::clicked, handler);
	currentTabLayout->addRow(button);
}

// View parameter updating methods

void Window::updateAngle(double dl, double db) {
	// dl and db are in degrees
	updatedViewParameters = true;
	config["l"] = static_cast<double>(config["l"]) + dl;
	l_field->setText(QString::number(static_cast<double>(config["l"])));
	config["b"] = static_cast<double>(config["b"]) + db;
	b_field->setText(QString::number(static_cast<double>(config["b"])));
}

void Window::updateViewSize(double factor) {
	// A positive factor means the view should become larger
	updatedViewParameters = true;
	config["view_width"] = static_cast<double>(config["view_width"])*factor;
	view_width_field->setText(QString::number(static_cast<double>(config["view_width"])));
	config["view_height"] = static_cast<double>(config["view_height"])*factor;
	view_height_field->setText(QString::number(static_cast<double>(config["view_height"])));
}

// Handlers

void Window::readAndPreprocessDataCube() {

	// Reads datacube from file and does as much pre-processing as possible without making assumptions about the regular grid size and beyond
	// The datacube file format is assumed to be the same as the one of the big dataset from the Frontiers article

	// Delete unnecessary objects
	if (FMO != NULL)
		delete FMO;
	if (renderer != NULL)
		delete renderer;

	// Initialization
	logProcessStart();
	int dims = static_cast<int>(config["dims"]);
	int nvars = static_cast<int>(config["nvars"]);
	int nx = static_cast<int>(config["nx"]);
	int ny = static_cast<int>(config["ny"]);
	int nz = static_cast<int>(config["nz"]);
	int ng = static_cast<int>(config["ng"]);
	FMO = new FoMo::FoMoObject(dims);

	// Load DataCube

	// Initialization
	std::string identity;
	std::ifstream in(static_cast<std::string const&>(config["datacube_file"]));
	in >> identity;
	std::vector<FoMo::tcoord> grid;
	FoMo::tcoord xvec(nx);
	FoMo::tcoord yvec(ny);
	FoMo::tcoord zvec(nz);
	grid.resize(dims);
	std::vector<FoMo::tphysvar> allvar;
	allvar.resize(nvars);
	for (int i = 0; i < dims; i++)
		grid[i].resize(ng);
	for (int i = 0; i < nvars; i++)
		allvar[i].resize(ng);
	// Reading coordinates
	for (int j = 0; j < nx; j++)
		in >> xvec.at(j);
	for (int j = 0; j < ny; j++)
		in >> yvec.at(j);
	double tmpvar;
	for (int j = 0; j < nz*2; j++) {
		in >> tmpvar;
		if (j < nz) zvec.at(j) = tmpvar;
	}
	// Writing coordinates
	for (int m = 0; m < nz; m++) {
		for (int l = 0; l < ny; l++) {
			for (int k = 0; k < nx; k++) {
				int j = m*ny*nx + l*nx + k;
				grid[0][j] = xvec[k];
				grid[1][j] = yvec[l];
				grid[2][j] = zvec[m];
			}
		}
	}
	// Read data
	//const int bufferSize = 1024;
	//char *buffer = new char[bufferSize];
	//in.read(buffer, bufferSize);
	//int charsRead = in.gcount();
	//int bufferPos = 0;
	for (int i = 0; i < nvars - 1; i++) {
		for (int j = 0; j < ng; j++) {
			in >> allvar[i][j];// = parseFloat(in, buffer, bufferSize, bufferPos, charsRead);
		}
	}
	//delete[] buffer;
	in.close();
	for(int i = 0; i < nvars; i++)
		std::cout << i << " allvar " << allvar[0][i] << std::endl;
	// Convert data
	for (int j = 0; j < ng; j++) {
		// Convert from g/cm^3 to /cm^3 (assuming that mean molecular mass is 1)
		allvar[0][j] /= 1e3*GSL_CONST_MKSA_MASS_PROTON;
		// Convert speeds to m/s (from km/s)
		allvar[2][j] *= 1e3;
		allvar[3][j] *= 1e3;
		allvar[4][j] = 0;
	}
	// Construct DataCube
	FMO->setdata(grid, allvar);
	logProcessFinished("loading DataCube");

	// Construct Goftcube
	FMO->constructGoftcube(static_cast<std::string const&>(config["chianti_file"]), static_cast<std::string const&>(config["abundance_file"]),
			static_cast<int>(config["lambda_pixel"]) == 1 ? FoMo::FoMoObservationType::Imaging : FoMo::FoMoObservationType::Spectroscopic);
	logProcessFinished("constructing GoftCube");

	// Construct renderer
	renderer = new FoMo::RegularGridRendererWrapper(FMO->readgoftcubepointer());
	logProcessFinished("initializing renderer");

	rendererState = WINDOW_RENDERER_STATE_HAS_DATACUBE;

}

void Window::constructRegularGrid() {
	if (rendererState < WINDOW_RENDERER_STATE_HAS_DATACUBE) {
		log("Please load a DataCube before attempting to construct a regular grid.");
	} else {
		// Only construct the regular grid if correct state has been reached
		logProcessStart();
		renderer->constructRegularGrid(static_cast<int>(config["gridx"]), static_cast<int>(config["gridy"]), static_cast<int>(config["gridz"]),
				static_cast<double>(config["max_distance_x"]), static_cast<double>(config["max_distance_y"]), static_cast<double>(config["max_distance_z"]));
		rendererState = WINDOW_RENDERER_STATE_HAS_REGULAR_GRID;
		logProcessFinished("constructing regular grid");
	}
}

void Window::updateRenderingSettings() {
	if (rendererState < WINDOW_RENDERER_STATE_HAS_REGULAR_GRID) {
		log("Please construct a regular grid before attempting to set the rendering settings.");
	} else {
		// Only set the rendering settings if correct state has been reached
		logProcessStart();
		int x_pixel = static_cast<int>(config["x_pixel"]);
		int y_pixel = static_cast<int>(config["y_pixel"]);
		// Update the renderer
		renderer->setRenderingSettings(x_pixel, y_pixel, static_cast<int>(config["lambda_pixel"]), static_cast<double>(config["lambda_width"]), FoMo::RegularGridRendererDisplayMode::IntegratedIntensity,
				static_cast<double>(config["max_intensity"]));
		// Update our own buffer and label sizes
		if (imageBuffer != NULL)
			delete[] imageBuffer;
		imageBuffer = new unsigned char[x_pixel*y_pixel];
		view->setFixedWidth(x_pixel);
		view->setFixedHeight(y_pixel);
		updatedRenderingSettings = true;
		rendererState = WINDOW_RENDERER_STATE_HAS_RENDERING_SETTINGS;
		logProcessFinished("updating rendering settings");
	}
}

void Window::renderToFile() {
	if (rendererState < WINDOW_RENDERER_STATE_HAS_RENDERING_SETTINGS) {
		log("Please set the rendering settings before attempting to render to a file");
	} else {
		// Only render to file if correct state has been reached
		logProcessStart();
		renderer->renderToCube(toRadians(static_cast<double>(config["l"])), toRadians(static_cast<double>(config["b"])), static_cast<double>(config["view_width"]), static_cast<double>(config["view_height"]),
				static_cast<std::string const&>(config["outfile"]));
		logProcessFinished("rendering to file");
	}
}

void Window::renderToView() {
	if (rendererState == WINDOW_RENDERER_STATE_HAS_RENDERING_SETTINGS) {
		// Only render to the view if correct state has been reached
		if (updatedRenderingSettings || updatedViewParameters) {
			logProcessStart();
			// Only re-render if any relevant settings have been updated
			int x_pixel = static_cast<int>(config["x_pixel"]);
			int y_pixel = static_cast<int>(config["y_pixel"]);
			renderer->renderToBuffer(toRadians(static_cast<double>(config["l"])), toRadians(static_cast<double>(config["b"])), static_cast<double>(config["view_width"]),
					static_cast<double>(config["view_height"]), imageBuffer);
			image = QImage(imageBuffer, x_pixel, y_pixel, 1*x_pixel, QImage::Format_Grayscale8);
			view->setPixmap(QPixmap::fromImage(image));
			updatedRenderingSettings = false;
			updatedViewParameters = false;
			logProcessFinished("rendering frame");
		}
	}
}

// Helper methods

inline float Window::parseFloat(std::ifstream &in, char *buffer, const int bufferSize, int &bufferPos, int &charsRead) {

	// Read one floating-point number from the stream
	// Manually implemented because this seems to be ~4 times faster than in >> allvar[i][j]
	// atof is already twice as fast but dependent on local settings (decimal point vs comma, ...)

	// buffer contains the last bufferSize characters before the current position in the input stream, or less if the end of the stream has been reached
	// bufferPos indicates our current position in the buffer
	// charsRead indicates the amount of characters read from the input stream into the buffer

	// Initialization
	const float powersOfTen[] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
	bool stop = false;
	bool foundDigits = false;
	int fractionalIndex = -1;
	bool inExponent = false;
	int mantissaSign = 1;
	int exponentSign = 1;
	float mantissa = 0;
	int exponent = 0;

	// Iterate through characters
	while(!stop) {
		if (bufferPos >= charsRead) {
			// Read more data from file
			in.read(buffer, bufferSize);
			charsRead = in.gcount();
			bufferPos = 0;
			if (charsRead == 0) {
				// File has ended
				break;
			}
		}
		char character = buffer[bufferPos++];
		// Handle every possible character
		switch (character) {
		case '0': case '1': case '2': case '3':
		case '4': case '5': case '6': case '7':
		case '8': case '9':
		{
			// Found digit, handle different cases
			foundDigits = true;
			int value = character - '0';
			if (!inExponent) {
				// In mantissa
				if (fractionalIndex == -1) {
					// In whole part
					mantissa = 10*mantissa + value;
				} else {
					// In fractional part
					mantissa += value*powersOfTen[fractionalIndex++];
				}
			} else {
				// In exponent
				exponent = 10*exponent + value;
			}
			break;
		}
		case '-':
			// Make current number negative
			if (!inExponent) {
				mantissaSign = -1;
			} else {
				exponentSign = -1;
			}
			break;
		case '.':
			// Entering fractional part
			fractionalIndex = 0;
			break;
		case 'e':
			// Entering exponent
			inExponent = 1;
			break;
		default:
			// If we reach any other character after finding digits, the number has ended
			stop = foundDigits;
		}
	}

	return mantissaSign*mantissa*pow(10, exponentSign*exponent);

}

inline float Window::toRadians(float degrees) {
	return degrees/180.0*M_PI;
}

inline std::chrono::time_point<std::chrono::high_resolution_clock> Window::time_now() {
	return std::chrono::high_resolution_clock::now();
}

inline void Window::log(const std::string &text) {
	logger->log(text);
}

inline void Window::logProcessStart() {
	start = time_now();
}

inline void Window::logProcessFinished(std::string name) {
	// Also resets start
	log(std::string("Finished ") + name + std::string(" in ")
		+ std::to_string(std::round(std::chrono::duration<double>(time_now() - start).count()*1000000)/1000000.0) + std::string(" seconds."));
	start = time_now();
}

void Window::printTest() {
	time += 0.016;
	//delete image;
	int width = 512;
	int height = 512;
	//image = new QImage(width, height, QImage::Format_ARGB32);
	for(int y = 0; y < height; y++) {
		for(int x = 0; x < width; x++) {
			int value = int((time*512 + x + y)/512.0*256)%256;
			image.setPixel(x, y, qRgb(value, value, value));
		}
	}
	//auto start = std::chrono::high_resolution_clock::now();
	view->setPixmap(QPixmap::fromImage(image));
	//auto duration = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	//log(std::to_string(duration));
}
