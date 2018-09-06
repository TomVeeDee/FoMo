#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <functional>

#include <math.h>

#include <QApplication>
#include <QStyle>
#include <QDesktopWidget>
#include <QGridLayout>
#include <QFormLayout>
#include <QTimer>
#include <QImage>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>

#include "window.h"
#include "ui_mainwindow.h"
#include "FoMo.h"

Window::Window(QApplication &app, QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow) {

	// Initializes this window
	// Caller is still responsible for calling show() when necessary

	// General
	ui->setupUi(this);
	resize(1280, 720);
	setWindowTitle("FoMo rendering GUI");
	setGeometry(QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(), app.desktop()->availableGeometry())); // Centers window on screen

	// Load config
	std::ifstream i(configPath);
	i >> config;
	i.close();

	// FoMo
	FMO = new FoMo::FoMoObject(3);

	// Top-level layout
	QWidget *centralWidget = new QWidget(this);
	QGridLayout *centralLayout = new QGridLayout(); // No need to assign parent here, setLayout already assigns ownership
	centralWidget->setLayout(centralLayout);
	setCentralWidget(centralWidget);

	// Image
	image = new QImage(512, 512, QImage::Format_ARGB32);
	for(int y = 0; y < 512; y++) {
		for(int x = 0; x < 512; x++) {
			int value = int(255*(time + x/512.0 + y/512.0))%256;
			image->setPixel(x, y, qRgb(value, value, value));
		}
	}
	label = new QLabel(centralWidget);
	label->setPixmap(QPixmap::fromImage(*image));
	centralLayout->addWidget(label, 0, 0, -1, 1, Qt::Alignment(Qt::AlignCenter));

	// Control panel initialization
	controlPanel = new QTabWidget(centralWidget);
	centralLayout->addWidget(controlPanel, 0, 1);

	// Datacube
	createTab("Datacube");
	addStringField("Datacube file", "datacube_file");
	addStringField("CHIANTI file", "chianti_file");
	addStringField("Abundance file", "abundance_file");
	addIntegerField("Dimensions", "dims", 2, 3);
	addIntegerField("Variables", "nvars", 1);
	addIntegerField("X points", "nx", 2);
	addIntegerField("Y points", "ny", 2);
	addIntegerField("Z points", "nz", 1);
	addIntegerField("Data points", "ng", 1);
	addButton("Load and preprocess datacube", [this] () {
		//TODO: Load datacube and preprocess (make GoftCube, build RTree?)
	});

	// Regular grid
	createTab("Regular grid");
	addIntegerField("X points", "gridx", 2);
	addIntegerField("Y points", "gridy", 2);
	addIntegerField("Z points", "gridz", 1);
	addButton("Construct regular grid", [this] () {
		//TODO: Construct regular grid
	});

	// Rendering
	createTab("Rendering");
	addIntegerField("X pixels", "x_pixel", 2);
	addIntegerField("Y pixels", "y_pixel", 2);
	addIntegerField("Lambda pixels", "lambda_pixel", 1);
	addDoubleField("View width (Mm)", "view_width");
	addDoubleField("View height (Mm)", "view_height");
	addDoubleField("Lambda width (m/s)", "lambda_width");
	addDoubleField("l (degrees)", "l");
	addDoubleField("b (degrees)", "b");
	addDoubleField("Frames/second", "fps");
	addButton("Update rendering settings", [this] () {
		//TODO: Update rendering settings
	});
	addStringField("Output file", "outfile");
	addButton("Render to file", [this] () {
		//TODO: Render to file
	});

	// Save config button
	QPushButton *save_config = new QPushButton("Save configuration to file", centralWidget);
	connect(save_config, &QPushButton::clicked, [this] () {
		std::ofstream o(this->configPath);
		o << std::setw(4) << this->config;
		o.close();
	});
	centralLayout->addWidget(save_config, 1, 1);

	// Log
	log = new Logger(centralWidget);
	centralLayout->addWidget(log, 2, 1);

	// Log + image test
	QTimer *timer = new QTimer(log);
	timer->start(33);
	connect(timer, &QTimer::timeout, this, &Window::PrintTest);

}

Window::~Window() {
	// Destructor, most destruction is handled by the Qt parent-child hierarchy
	delete FMO;
	delete image;
	delete ui;
}

// GUI creation methods

void Window::createTab(QString title) {
	// Creates a new tab with the given title
	currentTab = new QWidget(controlPanel);
	controlPanel->addTab(currentTab, title);
	currentTabLayout = new QFormLayout(); // No need to assign parent here, setLayout already assigns ownership
	currentTab->setLayout(currentTabLayout);
}

void Window::addStringField(QString title, std::string name) {
	// Creates a labelled string field connected to the config
	QLineEdit *lineEdit = new QLineEdit(QString::fromStdString(config[name]), currentTab);
	connect(lineEdit, &QLineEdit::textEdited, [this, name] (const QString &text) {
		this->config[name] = text.toStdString();
	});
	currentTabLayout->addRow(title, lineEdit);
}

void Window::addIntegerField(QString title, std::string name, int min, int max) {
	// max is inclusive
	QLineEdit *lineEdit = new QLineEdit(QString::number((int) config[name]), currentTab);
	lineEdit->setValidator(new QIntValidator(min, max, lineEdit));
	connect(lineEdit, &QLineEdit::textEdited, [this, name, lineEdit] (const QString &text) {
		if (lineEdit->hasAcceptableInput())
			this->config[name] = text.toInt();
	});
	currentTabLayout->addRow(title, lineEdit);
}

void Window::addDoubleField(QString title, std::string name, double min, double max) {
	QLineEdit *lineEdit = new QLineEdit(QString::number((double) config[name]), currentTab);
	lineEdit->setValidator(new QDoubleValidator(min, max, 10, lineEdit));
	connect(lineEdit, &QLineEdit::textEdited, [this, name, lineEdit] (const QString &text) {
		if (lineEdit->hasAcceptableInput())
			this->config[name] = text.toDouble();
	});
	currentTabLayout->addRow(title, lineEdit);
}

void Window::addButton(QString title, std::function<void (void)> handler) {
	QPushButton *button = new QPushButton(title, currentTab);
	connect(button, &QPushButton::clicked, handler);
	currentTabLayout->addRow(button);
}

void Window::PrintTest() {
	time += 0.033;
	delete image;
	int width = 512;
	int height = 512;
	image = new QImage(width, height, QImage::Format_ARGB32);
	for(int y = 0; y < height; y++) {
		for(int x = 0; x < width; x++) {
			int value = int((time*512 + x + y)/512.0*256)%256;
			image->setPixel(x, y, qRgb(value, value, value));
		}
	}
	auto start = std::chrono::high_resolution_clock::now();
	label->setPixmap(QPixmap::fromImage(*image));
	auto duration = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	log->log(std::to_string(duration));
}
