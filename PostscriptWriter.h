#pragma once

#include <fstream>

class Vector2;

/** @brief PostScriptファイルの書き出しクラス
 *  @attention
 *  comment(), show()で'%'や改行のエスケープを行っていないので注意
 */
class PostscriptWriter
{
public:
	explicit PostscriptWriter(const std::string& fileName);
	~PostscriptWriter();

	PostscriptWriter& moveTo(const Vector2& p);

	PostscriptWriter& lineTo(const Vector2& p);
	
	PostscriptWriter& stroke();

	PostscriptWriter& fill();

	PostscriptWriter& setFontBatch(char name[], unsigned int size);

	PostscriptWriter& show(const char s[]);

	PostscriptWriter& comment(const char s[]);

	PostscriptWriter& showpage();

private:
	void close();

	std::ostream& header();

	std::ostream& point(const Vector2& p);

private:
	std::ofstream fout_;
};
