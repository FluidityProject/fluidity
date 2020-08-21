#ifndef QVISH5PARTRANDOMSAMPLEWINDOW_H
#define QVISH5PARTRANDOMSAMPLEWINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>

class H5PartRandomSampleAttributes;
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
// Class: QvisH5PartRandomSampleWindow
//
// Purpose:
//   Defines QvisH5PartRandomSampleWindow class.
//
// Notes:      This class was automatically generated!

// Programmer: xml2window
// Creation:   Thu Mar 16 10:26:55 PDT 2006
//
// Modifications:
//
// ****************************************************************************

class QvisH5PartRandomSampleWindow : public QvisOperatorWindow
{
	Q_OBJECT
public:
	QvisH5PartRandomSampleWindow(const int type,
	                             H5PartRandomSampleAttributes *subj,
	                             const char *caption = 0,
	                             const char *shortName = 0,
	                             QvisNotepadArea *notepad = 0);
	virtual ~QvisH5PartRandomSampleWindow();
	virtual void CreateWindowContents();
protected:
	void UpdateWindow(bool doAll);
	virtual void GetCurrentValues(int which_widget);
private slots:
	void factorProcessText();
private:
	QLineEdit *factor;
	QLabel *factorLabel;

	H5PartRandomSampleAttributes *atts;
};



#endif
