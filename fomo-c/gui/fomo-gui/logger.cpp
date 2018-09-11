#include <QPlainTextEdit>
#include <QScrollBar>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <iostream>

#include "logger.h"

Logger::Logger(QWidget *parent) : QPlainTextEdit(parent) {
	setReadOnly(true);
	empty = true;
}

void Logger::log(const std::string &text) {
	std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	std::stringstream ss;
	ss << (empty ? "" : "\n") << "[" << std::put_time(std::localtime(&now), "%T") << "] " << text;
	moveCursor(QTextCursor::End);
	insertPlainText(QString::fromStdString(ss.str())); // Adds the text to the log without newline
	moveCursor(QTextCursor::End);
	empty = false;
}
