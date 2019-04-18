#ifndef QVISBINNERWINDOW_H
#define QVISBINNERWINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>

class Binner;
class QLabel;
class QCheckBox;
class QLineEdit;
class QSpinBox;
class QVBox;
class QButtonGroup;
class QvisColorTableButton;
class QvisOpacitySlider;
class QvisColorButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class QvisVariableButton;

// ****************************************************************************
// Class: QvisBinnerWindow
//
// Purpose:
//   Defines QvisBinnerWindow class.
//
// Notes:      This class was automatically generated!

// Programmer: xml2window
// Creation:   Thu Mar 30 12:05:26 PDT 2006
//
// Modifications:
//
// ****************************************************************************

class QvisBinnerWindow : public QvisOperatorWindow
{
	Q_OBJECT
public:
	QvisBinnerWindow(const int type,
	                 Binner *subj,
	                 const char *caption = 0,
	                 const char *shortName = 0,
	                 QvisNotepadArea *notepad = 0);
	virtual ~QvisBinnerWindow();
	virtual void CreateWindowContents();
protected:
	void UpdateWindow(bool doAll);
	virtual void GetCurrentValues(int which_widget);
private slots:
	void dim1ProcessText();
	void dim2ProcessText();
	void dim3ProcessText();
private:
	QLineEdit *dim1;
	QLineEdit *dim2;
	QLineEdit *dim3;
	QLabel *dim1Label;
	QLabel *dim2Label;
	QLabel *dim3Label;

	Binner *atts;
};



#endif
