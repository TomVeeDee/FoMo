#pragma once

#include <QPlainTextEdit>

class Logger : public QPlainTextEdit {

public:
	explicit Logger(QWidget *parent = 0);
	void log(const std::string &text);

private:
	bool empty = true;

};
