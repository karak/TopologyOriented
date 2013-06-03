#include "stdafx.h"
#include "PostscriptWriter.h"
#include "Vector2.h"


PostscriptWriter::PostscriptWriter(const std::string& fileName) : fout_(fileName.c_str())
{
	header();
}

PostscriptWriter::~PostscriptWriter()
{
	close();
}

/** -> lineto */
PostscriptWriter& PostscriptWriter::moveTo(const Vector2& p)
{
	point(p) << " moveto\n";
	return *this;
}

PostscriptWriter& PostscriptWriter::setFontBatch(char name[], unsigned int size)
{
	fout_
		<< '/' << name << " findfont\n"
		<< size << " scalefont\n"
		<< "setfont\n";
	return *this;
}

/** -> lineto | stroke */
PostscriptWriter& PostscriptWriter::lineTo(const Vector2& p)
{
	point(p) << " lineto\n";
	return *this;
}

PostscriptWriter& PostscriptWriter::stroke()
{
	fout_ << "stroke\n";
	return *this;
}

PostscriptWriter& PostscriptWriter::fill()
{
	fout_ << "fill\n";
	return *this;
}

PostscriptWriter& PostscriptWriter::showpage()
{
	fout_ << "showpage\n";
	return *this;
}
PostscriptWriter& PostscriptWriter::show(const char s[])
{
	fout_ << '(' << s << ')' << " show\n";
	return *this;
}

PostscriptWriter& PostscriptWriter::comment(const char s[])
{
	fout_ << '\n';
	bool isHeadOfLine = true;
	for (const char* p = s; *p != '\0'; ++p)
	{
		if (isHeadOfLine)  fout_ << '%';
		fout_ << *p;
		isHeadOfLine = (*p == '\n');
	}
	if (!isHeadOfLine)
		fout_ << '\n';

	return *this;
}

void PostscriptWriter::close()
{
	//ÅŒã‚É•K‚¸showpage‚Â‚¯‚é
	fout_ << "\n";
	showpage();
	fout_.close();
}

std::ostream& PostscriptWriter::header()
{
	fout_ << "%!\n";
	return fout_;
}

std::ostream& PostscriptWriter::point(const Vector2& p)
{
	return fout_ << p.x << ' ' << p.y;
}